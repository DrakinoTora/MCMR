#include "simulation_mg.hpp"
#include "physics.hpp"
#include "exporter.hpp"
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <pybind11/pybind11.h>
namespace py = pybind11;

SimulationMG::SimulationMG(int N, double x_world, double y_world,
                        const std::vector<double>& x_grid,
                        const std::vector<double>& y_grid,
                        const std::vector<std::vector<std::string>>& material_matrix,
                        const std::vector<std::vector<double>>& sources,
                        const std::vector<double>& circle_cx,
                        const std::vector<double>& circle_cy,
                        const std::vector<double>& circle_r,
                        const std::vector<std::string>& circle_material,
                        const std::vector<double>& circle_source,
                        int max_save,
                        const std::string& bc_top,
                        const std::string& bc_bot,
                        const std::string& bc_left,
                        const std::string& bc_right)
    : N_particles(N),
      grid(x_world, y_world, x_grid, y_grid, material_matrix,
           bc_left, bc_right, bc_bot, bc_top),
      max_history_save(max_save) {

    // -- geometry / source setup: identical logic to Simulation's constructor --
    int nx = grid.nx();
    int ny = grid.ny();
    n_grid_cells = nx * ny;

    if (static_cast<int>(sources.size()) != ny)
        throw std::invalid_argument(
            "sources row count (" + std::to_string(sources.size()) +
            ") must be same as len(y_grid)+1 = " + std::to_string(ny));

    int n_circles = static_cast<int>(circle_cx.size());
    if (static_cast<int>(circle_cy.size()) != n_circles ||
        static_cast<int>(circle_r.size()) != n_circles ||
        static_cast<int>(circle_material.size()) != n_circles ||
        static_cast<int>(circle_source.size()) != n_circles)
        throw std::invalid_argument(
            "circle_cx, circle_cy, circle_r, circle_material, circle_source must all have the same length");

    flat_source_weights.resize(n_grid_cells + n_circles);
    double total = 0.0;
    for (int row = 0; row < ny; ++row) {
        if (static_cast<int>(sources[row].size()) != nx)
            throw std::invalid_argument(
                "sources row " + std::to_string(row) +
                " col count must be same as len(x_grid)+1 = " + std::to_string(nx));

        int iy = ny - 1 - row;
        for (int col = 0; col < nx; ++col) {
            int ix = col;
            double w = sources[row][col];
            if (w < 0.0)
                throw std::invalid_argument("sources can't be negative value");
            flat_source_weights[ix * ny + iy] = w;
            total += w;
        }
    }

    for (int i = 0; i < n_circles; ++i) {
        if (circle_source[i] < 0.0)
            throw std::invalid_argument("circle source can't be negative value");
        grid.add_circle(circle_cx[i], circle_cy[i], circle_r[i], circle_material[i]);
        flat_source_weights[n_grid_cells + i] = circle_source[i];
        total += circle_source[i];
    }

    if (total <= 0.0)
        throw std::invalid_argument("sum of all sources (grid cells + circles) must be greater than 0");
}

void SimulationMG::set_group_data(
    int n_groups_,
    const std::map<std::string, std::vector<double>>& sigma_t_,
    const std::map<std::string, std::vector<double>>& sigma_s_,
    const std::map<std::string, std::vector<double>>& sigma_f_,
    const std::map<std::string, std::vector<double>>& nu_
) {
    if (n_groups_ <= 0)
        throw std::invalid_argument("n_groups must be positive");

    // sigma_f and nu are optional per material: a material missing from those maps
    // gets zeros (no fission), so run() can always index them by material symbol
    std::map<std::string, std::vector<double>> sf_full, nu_full;

    for (const auto& kv : sigma_t_) {
        const std::string& name = kv.first;
        if (static_cast<int>(kv.second.size()) != n_groups_)
            throw std::invalid_argument("sigma_t for material '" + name +
                                         "' must have exactly n_groups elements");
        if (!sigma_s_.count(name) || static_cast<int>(sigma_s_.at(name).size()) != n_groups_)
            throw std::invalid_argument("sigma_s missing or wrong length for material '" + name + "'");

        const auto& st = kv.second;
        const auto& ss = sigma_s_.at(name);

        std::vector<double> sf(n_groups_, 0.0);
        if (sigma_f_.count(name)) {
            if (static_cast<int>(sigma_f_.at(name).size()) != n_groups_)
                throw std::invalid_argument("sigma_f wrong length for material '" + name + "'");
            sf = sigma_f_.at(name);
        }

        std::vector<double> nu_m(n_groups_, 0.0);
        if (nu_.count(name)) {
            if (static_cast<int>(nu_.at(name).size()) != n_groups_)
                throw std::invalid_argument("nu wrong length for material '" + name + "'");
            nu_m = nu_.at(name);
        }

        for (int g = 0; g < n_groups_; ++g) {
            if (ss[g] + sf[g] > st[g] + 1e-12)
                throw std::invalid_argument("sigma_s + sigma_f exceeds sigma_t at group " + std::to_string(g) +
                                             " for material '" + name + "'");
            if (nu_m[g] < 0.0)
                throw std::invalid_argument("nu can't be negative at group " + std::to_string(g) +
                                             " for material '" + name + "'");
        }

        sf_full[name] = sf;
        nu_full[name] = nu_m;
    }

    n_groups = n_groups_;
    sigma_t = sigma_t_;
    sigma_s = sigma_s_;
    sigma_f = sf_full;
    nu = nu_full;
}

void SimulationMG::transport_one(double x, double y, double mu_x, double mu_y, int g, int ix, int iy,
                                 std::mt19937& gen,
                                 std::vector<double>* h_x, std::vector<double>* h_y) {
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
    const bool save_history = (h_x != nullptr && h_y != nullptr);

    bool alive = true;
    // Guards against a medium with zero absorption/fission probability inside
    // an all-reflective enclosure, where a neutron can never leak or die and
    // this loop would otherwise never return (and, since it never reaches a
    // Python-visible frame, Ctrl-C could not interrupt it either).
    constexpr long long kMaxCollisions = 5'000'000;
    long long n_collisions = 0;
    while (alive) {
        if ((++n_collisions % 100'000) == 0) {
            if (PyErr_CheckSignals() != 0) throw py::error_already_set();
            if (n_collisions >= kMaxCollisions)
                throw std::runtime_error(
                    "transport_one: exceeded " + std::to_string(kMaxCollisions) +
                    " collisions for a single neutron. This usually means the medium has "
                    "zero absorption/fission probability in some group while all boundaries "
                    "are reflective, so the neutron can never leak or die. Check sigma_a "
                    "(sigma_t - sigma_s - sigma_f) and the geometry's boundary conditions.");
        }
        const MaterialInfo& cur_mat = grid.material_at(x, y, ix, iy);
        const auto& st = sigma_t.at(cur_mat.symbol);
        const auto& ss = sigma_s.at(cur_mat.symbol);
        const auto& sf = sigma_f.at(cur_mat.symbol);

        double Sigma_t = st[g];
        double Sigma_s = ss[g];
        double Sigma_f = sf[g];
        double Sigma_a = Sigma_t - Sigma_s - Sigma_f; // pure capture

        double R = dist_R(gen);
        double d_coll = -std::log(R) / Sigma_t;

        Side hit_side;
        double d_surf = grid.distance_to_boundary(x, y, mu_x, mu_y, ix, iy, hit_side);

        double d = std::min(d_coll, d_surf);
        x += d * mu_x;
        y += d * mu_y;

        if (d_coll >= d_surf) {
            if (hit_side != Side::None) {
                BoundaryType bc = grid.bc_for_side(hit_side);
                if (bc == BoundaryType::Vacuum) {
                    results.transmission++;
                    alive = false;
                    break;
                } else {
                    if (hit_side == Side::Left || hit_side == Side::Right) mu_x = -mu_x;
                    else mu_y = -mu_y;
                    if (save_history) { h_x->push_back(x); h_y->push_back(y); }
                    continue;
                }
            }
            x += 1e-6 * mu_x;
            y += 1e-6 * mu_y;
            grid.find_index(x, y, ix, iy);
            continue;
        }

        if (save_history) { h_x->push_back(x); h_y->push_back(y); }

        // pick event by slicing [0, Sigma_t) into scatter / capture / fission
        double P = dist_R(gen) * Sigma_t;
        if (P < Sigma_s) {
            // isotropic scatter, down-scatter only: new group uniform in [0, g]
            double phi = dist_phi(gen);
            mu_x = std::cos(phi);
            mu_y = std::sin(phi);
            std::uniform_int_distribution<int> dist_newg(0, g);
            g = dist_newg(gen);
        } else if (P < Sigma_s + Sigma_a) {
            alive = false;
            results.absorp_by_material[cur_mat.symbol]++;
        } else {
            // fission: the parent dies here, its children go into the fission bank.
            // position = parent's position, group = parent's group,
            // direction = a fresh isotropic direction for EACH child
            alive = false;
            results.fission_by_material[cur_mat.symbol]++;

            int n_children = sample_fission_neutrons(nu.at(cur_mat.symbol)[g]);
            for (int k = 0; k < n_children; ++k) {
                double child_phi = dist_phi(gen);
                fission_bank.push_back({x, y, std::cos(child_phi), std::sin(child_phi), g});
            }
        }
    }
}

void SimulationMG::run() {
    if (n_groups <= 0)
        throw std::invalid_argument("set_group_data() must be called before run()");

    auto start_time = std::chrono::high_resolution_clock::now();

    // Reset state left over from a PREVIOUS run() call on this same object --
    // both the Tally (so results don't silently accumulate across two runs)
    // and the fission bank (an interrupted earlier run() must not leak
    // neutrons into this one; fission_bank.clear() below already covered
    // that half, results = Tally{} covers the rest).
    results = Tally{};
    fission_bank.clear();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
    std::discrete_distribution<int> region_picker(flat_source_weights.begin(), flat_source_weights.end());
    std::uniform_int_distribution<int> dist_g0(0, n_groups - 1);

    py::print("====================================================");
    py::print(" MCMR Simulation Engine (multi-group)");
    py::print(" Total Particle   :", N_particles);
    py::print(" Number of Groups :", n_groups);
    py::print("====================================================");
    py::print("running...");
    py::module_::import("sys").attr("stdout").attr("flush")();

    const int ny = grid.ny();

    // Progress line: [source particle in progress / N] + current fission bank size.
    // Throttled by time (not by particle count) because a supercritical run can sit on ONE
    // source particle forever -- the bank size is what tells the user that is happening.
    auto last_print = std::chrono::steady_clock::now();
    auto print_progress = [&](int current, bool force) {
        auto now = std::chrono::steady_clock::now();
        if (!force && now - last_print < std::chrono::milliseconds(100)) return;
        last_print = now;

        // lets Ctrl+C / "Interrupt kernel" stop a run that never ends
        if (PyErr_CheckSignals() != 0) throw py::error_already_set();

        int percent = N_particles > 0
            ? static_cast<int>(static_cast<long long>(current) * 100 / N_particles)
            : 100;

        std::stringstream ss;
        ss << "\rProgress: [" << current << "/" << N_particles << "] (" << percent << "%)"
           << " | Fission bank: " << fission_bank.size() << "      ";  // trailing spaces wipe a longer old line

        py::print(ss.str(), py::arg("end") = "");
        py::module_::import("sys").attr("stdout").attr("flush")();
    };

    for (int p = 0; p < N_particles; ++p) {
        print_progress(p + 1, p == 0);

        bool save_history = p < max_history_save;
        std::vector<double> h_x, h_y;

        int flat_idx = region_picker(gen);
        double x, y;
        int ix, iy;

        if (flat_idx < n_grid_cells) {
            ix = flat_idx / ny;
            iy = flat_idx % ny;
            const Region& r = grid.region_at(ix, iy);
            x = r.x1 + dist_R(gen) * (r.x2 - r.x1);
            y = r.y1 + dist_R(gen) * (r.y2 - r.y1);
        } else {
            const CircleRegion& c = grid.circle_at_index(flat_idx - n_grid_cells);
            double u = dist_R(gen);
            double theta = dist_phi(gen);
            double rr = c.r * std::sqrt(u); // uniform over the disk area, see region-sampling note
            x = c.cx + rr * std::cos(theta);
            y = c.cy + rr * std::sin(theta);
            grid.find_index(x, y, ix, iy);
        }

        double phi = dist_phi(gen);
        double mu_x = std::cos(phi);
        double mu_y = std::sin(phi);
        int g = dist_g0(gen); // group born uniformly

        if (save_history) { h_x.push_back(x); h_y.push_back(y); }

        transport_one(x, y, mu_x, mu_y, g, ix, iy, gen,
                      save_history ? &h_x : nullptr,
                      save_history ? &h_y : nullptr);

        // trajectories are saved for source particles only (max_history_save keeps its meaning)
        if (save_history) {
            results.x_history.push_back(h_x);
            results.y_history.push_back(h_y);
        }

        // Fission bank: do NOT move on to the next source particle until the bank is empty.
        // Banked neutrons are transported exactly like any other neutron; a fission in here
        // just appends more members to this same queue.
        while (!fission_bank.empty()) {
            print_progress(p + 1, false);

            BankedNeutron n = fission_bank.front();  // copy first: transport_one() pushes to the bank
            fission_bank.pop_front();

            int bix, biy;
            grid.find_index(n.x, n.y, bix, biy);
            transport_one(n.x, n.y, n.mu_x, n.mu_y, n.g, bix, biy, gen, nullptr, nullptr);
        }
    }

    print_progress(N_particles, true);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    results.time_taken = diff.count();

    py::print("\n====================================================");
    std::stringstream ss_end;
    ss_end << " Simulation End in " << std::fixed << std::setprecision(4) << results.time_taken << " seconds.";
    py::print(ss_end.str());
    py::print("====================================================\n");
    py::module_::import("sys").attr("stdout").attr("flush")();

    // run() no longer writes to disk on its own -- call sim.export_xml("your_file.xml")
    // explicitly afterward if/when you actually want a file. See export_xml() below.
}

void SimulationMG::export_xml(const std::string& filename) {
    export_to_xml(results, filename);
}
