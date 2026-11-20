#include "simulation_mg.hpp"
#include "physics.hpp"
#include "rng.hpp"
#include "exporter.hpp"
#include <algorithm>
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
                        int tally_res_x_,
                        int tally_res_y_,
                        const std::string& bc_top,
                        const std::string& bc_bot,
                        const std::string& bc_left,
                        const std::string& bc_right)
    : N_particles(N),
      grid(x_world, y_world, x_grid, y_grid, material_matrix,
           bc_left, bc_right, bc_bot, bc_top),
      max_history_save(max_save),
      // 0 (unset) falls back to one tally cell per world unit -- e.g. a
      // 60x50 world defaults to a 60x50 tally grid (1x1 cells)
      tally_res_x(tally_res_x_ > 0 ? tally_res_x_ : std::max(1, static_cast<int>(std::lround(x_world)))),
      tally_res_y(tally_res_y_ > 0 ? tally_res_y_ : std::max(1, static_cast<int>(std::lround(y_world)))) {

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
    const std::map<std::string, std::vector<double>>& nu_,
    const std::map<std::string, std::vector<double>>& chi_fission_,
    const std::map<std::string, std::vector<std::vector<double>>>& chi_scatter_
) {
    if (n_groups_ <= 0)
        throw std::invalid_argument("n_groups must be positive");

    // sigma_f and nu are optional per material: a material missing from those maps
    // gets zeros (no fission), so run() can always index them by material symbol
    std::map<std::string, std::vector<double>> sf_full, nu_full;
    // chi_fission is optional too, but "absent" means something different from
    // zeros (inherit the parent's group), so only materials that really give
    // a chi get an entry here
    std::map<std::string, std::vector<double>> chi_cdf_full;
    // same for chi_scatter: only materials that really give a matrix get an entry
    // ("absent" = the older uniform [0, g] rule), stored as one cdf per row
    std::map<std::string, std::vector<std::vector<double>>> chi_scatter_cdf_full;

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

        auto chi_it = chi_fission_.find(name);
        if (chi_it != chi_fission_.end() && !chi_it->second.empty()) {
            const auto& chi = chi_it->second;
            if (static_cast<int>(chi.size()) != n_groups_)
                throw std::invalid_argument("chi_fission for material '" + name +
                                             "' must have exactly n_groups elements");
            double total = 0.0;
            for (int g = 0; g < n_groups_; ++g) {
                if (!std::isfinite(chi[g]) || chi[g] < 0.0)
                    throw std::invalid_argument("chi_fission must be finite and non-negative, got " +
                                                 std::to_string(chi[g]) + " at group " + std::to_string(g) +
                                                 " for material '" + name + "'");
                total += chi[g];
            }
            if (!(total > 0.0))
                throw std::invalid_argument("chi_fission for material '" + name +
                                             "' must have a positive sum (can't normalize all zeros)");

            std::vector<double> cdf(n_groups_);
            double running = 0.0;
            for (int g = 0; g < n_groups_; ++g) {
                running += chi[g];
                cdf[g] = running / total;
            }
            cdf.back() = 1.0;  // kill any rounding gap, so a draw u < 1 always lands on a group
            chi_cdf_full[name] = std::move(cdf);
        }

        auto cs_it = chi_scatter_.find(name);
        if (cs_it != chi_scatter_.end() && !cs_it->second.empty()) {
            const auto& mat = cs_it->second;
            if (static_cast<int>(mat.size()) != n_groups_)
                throw std::invalid_argument("chi_scatter for material '" + name +
                                             "' must have exactly n_groups rows (one per group before the scatter)");

            std::vector<std::vector<double>> cdf_rows(n_groups_);
            for (int g = 0; g < n_groups_; ++g) {
                const auto& row = mat[g];
                if (static_cast<int>(row.size()) != n_groups_)
                    throw std::invalid_argument("chi_scatter row " + std::to_string(g) + " for material '" + name +
                                                 "' must have exactly n_groups elements");
                double total = 0.0;
                for (int j = 0; j < n_groups_; ++j) {
                    if (!std::isfinite(row[j]) || row[j] < 0.0)
                        throw std::invalid_argument("chi_scatter must be finite and non-negative, got " +
                                                     std::to_string(row[j]) + " at row " + std::to_string(g) +
                                                     ", column " + std::to_string(j) + " for material '" + name + "'");
                    total += row[j];
                }

                std::vector<double> cdf(n_groups_);
                if (total > 0.0) {
                    double running = 0.0;
                    for (int j = 0; j < n_groups_; ++j) {
                        running += row[j];
                        cdf[j] = running / total;
                    }
                    cdf.back() = 1.0;  // kill any rounding gap, so a draw u < 1 always lands on a group
                } else {
                    // nothing to normalize. Only fine if neutrons never scatter out of group g
                    if (ss[g] > 0.0)
                        throw std::invalid_argument("chi_scatter row " + std::to_string(g) + " for material '" + name +
                                                     "' sums to zero, but sigma_s[" + std::to_string(g) +
                                                     "] > 0 (a scattered neutron would have nowhere to go)");
                    // never sampled; park a harmless "stay in group g" row
                    for (int j = 0; j < n_groups_; ++j) cdf[j] = (j >= g) ? 1.0 : 0.0;
                }
                cdf_rows[g] = std::move(cdf);
            }
            chi_scatter_cdf_full[name] = std::move(cdf_rows);
        }
    }

    n_groups = n_groups_;
    sigma_t = sigma_t_;
    sigma_s = sigma_s_;
    sigma_f = sf_full;
    nu = nu_full;
    chi_cdf = chi_cdf_full;
    chi_scatter_cdf = chi_scatter_cdf_full;
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
        double x_prev = x, y_prev = y;
        x += d * mu_x;
        y += d * mu_y;
        // Sigma_f / Sigma_a of THIS segment's material and group (a segment never
        // spans two materials); Sigma_a clamped so inconsistent data can't score negative
        results.add_flux_track(x_prev, y_prev, x, y, g, Sigma_f, std::max(0.0, Sigma_a));

        if (d_coll >= d_surf) {
            if (hit_side != Side::None) {
                BoundaryType bc = grid.bc_for_side(hit_side);
                if (bc == BoundaryType::Vacuum) {
                    switch (hit_side) {
                        case Side::Left:   results.leak_left++;   break;
                        case Side::Right:  results.leak_right++;  break;
                        case Side::Bottom: results.leak_bottom++; break;
                        case Side::Top:    results.leak_top++;    break;
                        default: break;
                    }
                    results.G_leak.push_back(g);
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
            // isotropic scatter. New group: row g of this material's chi_scatter if it has
            // one (a draw from that row's cdf -- any group is allowed), otherwise the
            // older rule, uniform in [0, g]
            double phi = dist_phi(gen);
            mu_x = std::cos(phi);
            mu_y = std::sin(phi);

            const auto cs_it = chi_scatter_cdf.find(cur_mat.symbol);
            if (cs_it != chi_scatter_cdf.end()) {
                const std::vector<double>& row_cdf = cs_it->second[g];
                // first group whose cumulative probability exceeds u. Columns with
                // probability 0 repeat the previous cdf value, so they are never picked.
                double u = dist_R(gen);
                g = static_cast<int>(std::upper_bound(row_cdf.begin(), row_cdf.end(), u) - row_cdf.begin());
            } else {
                std::uniform_int_distribution<int> dist_newg(0, g);
                g = dist_newg(gen);
            }
        } else if (P < Sigma_s + Sigma_a) {
            alive = false;
            results.add_absorp_hit(x, y);
        } else {
            // fission: the parent dies here, its children go into the fission bank.
            // position = parent's position,
            // direction = a fresh isotropic direction for EACH child,
            // group = drawn for EACH child from this material's chi_fission
            //         (independent of the parent's group g); if the material
            //         has no chi_fission, the child just inherits g
            alive = false;
            results.add_fission_hit(x, y);

            const auto chi_it = chi_cdf.find(cur_mat.symbol);
            const std::vector<double>* cdf = (chi_it != chi_cdf.end()) ? &chi_it->second : nullptr;

            int n_children = sample_fission_neutrons(nu.at(cur_mat.symbol)[g] / k_current, gen);
            for (int k = 0; k < n_children; ++k) {
                double child_phi = dist_phi(gen);
                int child_g = g;
                if (cdf) {
                    // first group whose cumulative probability exceeds u. Groups with
                    // chi == 0 repeat the previous cdf value, so they are never picked.
                    double u = dist_R(gen);
                    child_g = static_cast<int>(std::upper_bound(cdf->begin(), cdf->end(), u) - cdf->begin());
                }
                fission_bank.push_back({x, y, std::cos(child_phi), std::sin(child_phi), child_g});
            }
        }
    }
}

void SimulationMG::run(int n_gen, double k, std::optional<std::uint64_t> seed) {
    if (n_groups <= 0)
        throw std::invalid_argument("set_group_data() must be called before run()");
    if (n_gen <= 0)
        throw std::invalid_argument("n_gen must be positive");
    if (!(k > 0.0))
        throw std::invalid_argument("k must be positive");

    auto start_time = std::chrono::high_resolution_clock::now();

    // Reset state left over from a PREVIOUS run() call on this same object --
    // both the Tally (so results don't silently accumulate across two runs)
    // and the fission bank (an interrupted earlier run() must not leak
    // neutrons into this one; fission_bank.clear() below already covered
    // that half, results = Tally{} covers the rest).
    results = Tally{};
    results.seed = resolve_seed(seed);
    results.init_region_tally(tally_res_x, tally_res_y, grid.world_max_x(), grid.world_max_y(), n_groups);
    fission_bank.clear();
    k_current = k;

    std::mt19937 gen = make_rng(results.seed);
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
    std::discrete_distribution<int> region_picker(flat_source_weights.begin(), flat_source_weights.end());
    std::uniform_int_distribution<int> dist_g0(0, n_groups - 1);

    py::print("====================================================");
    py::print(" MCMR Simulation Engine (multi-group)");
    py::print(" Particles / Generation :", N_particles);
    py::print(" Generations            :", n_gen);
    py::print(" Number of Groups       :", n_groups);
    py::print(" Seed                   :", results.seed);
    py::print("====================================================");
    py::print("running...");
    py::module_::import("sys").attr("stdout").attr("flush")();

    const int ny = grid.ny();

    // Progress line: current generation + source particle in progress / N + current
    // fission bank size (children accumulated so far THIS generation). Throttled by
    // time (not by particle count) because a supercritical run can sit on one
    // neutron's collision chain for a while -- the bank size is what tells the
    // user something is still happening.
    auto last_print = std::chrono::steady_clock::now();
    auto print_progress = [&](int current, int gen_idx, bool force) {
        auto now = std::chrono::steady_clock::now();
        if (!force && now - last_print < std::chrono::milliseconds(100)) return;
        last_print = now;

        // lets Ctrl+C / "Interrupt kernel" stop a run that never ends
        if (PyErr_CheckSignals() != 0) throw py::error_already_set();

        int percent = N_particles > 0
            ? static_cast<int>(static_cast<long long>(current) * 100 / N_particles)
            : 100;

        std::stringstream ss;
        ss << "\rGen " << gen_idx << "/" << n_gen
           << " | Particle: [" << current << "/" << N_particles << "] (" << percent << "%)"
           << " | Fission bank: " << fission_bank.size() << "      ";  // trailing spaces wipe a longer old line

        py::print(ss.str(), py::arg("end") = "");
        py::module_::import("sys").attr("stdout").attr("flush")();
    };

    // this generation's population, transported one by one below. Generation 1:
    // sampled from the spatial/group source distribution. Generation >1: N_particles
    // draws WITH REPLACEMENT from the previous generation's fission_bank (a fixed
    // population size regardless of whether the system is sub/super/critical).
    std::vector<BankedNeutron> population;
    population.reserve(N_particles);

    int gen_idx = 1;
    for (; gen_idx <= n_gen; ++gen_idx) {
        population.clear();

        if (gen_idx == 1) {
            for (int i = 0; i < N_particles; ++i) {
                int flat_idx = region_picker(gen);
                double x, y;

                if (flat_idx < n_grid_cells) {
                    int gix = flat_idx / ny;
                    int giy = flat_idx % ny;
                    const Region& r = grid.region_at(gix, giy);
                    x = r.x1 + dist_R(gen) * (r.x2 - r.x1);
                    y = r.y1 + dist_R(gen) * (r.y2 - r.y1);
                } else {
                    const CircleRegion& c = grid.circle_at_index(flat_idx - n_grid_cells);
                    double u = dist_R(gen);
                    double theta = dist_phi(gen);
                    double rr = c.r * std::sqrt(u); // uniform over the disk area
                    x = c.cx + rr * std::cos(theta);
                    y = c.cy + rr * std::sin(theta);
                }

                double phi = dist_phi(gen);
                int g = dist_g0(gen); // group born uniformly
                population.push_back({x, y, std::cos(phi), std::sin(phi), g});
            }
        } else {
            if (fission_bank.empty()) {
                py::print("\nChain died out before generation", gen_idx,
                           "(no fissions in the previous generation) -- stopping early.");
                py::module_::import("sys").attr("stdout").attr("flush")();
                break;
            }
            std::uniform_int_distribution<size_t> dist_bank(0, fission_bank.size() - 1);
            for (int i = 0; i < N_particles; ++i) {
                const BankedNeutron& src = fission_bank[dist_bank(gen)]; // WITH replacement
                double phi = dist_phi(gen); // fresh isotropic direction for each draw
                population.push_back({src.x, src.y, std::cos(phi), std::sin(phi), src.g});
            }
        }

        // this generation's bank has now been fully folded into `population` (or,
        // for generation 1, was already empty) -- clear it so it only accumulates
        // THIS generation's fresh fission events from here on
        fission_bank.clear();

        // the k this generation runs with (recorded in the XML, not printed)
        results.k_used.push_back(k_current);

        for (int p = 0; p < N_particles; ++p) {
            print_progress(p + 1, gen_idx, p == 0);

            bool save_history = p < max_history_save;
            std::vector<double> h_x, h_y;

            const BankedNeutron& sp = population[p];
            int ix, iy;
            grid.find_index(sp.x, sp.y, ix, iy);

            results.G_born.push_back(sp.g);
            if (save_history) { h_x.push_back(sp.x); h_y.push_back(sp.y); }

            transport_one(sp.x, sp.y, sp.mu_x, sp.mu_y, sp.g, ix, iy, gen,
                          save_history ? &h_x : nullptr,
                          save_history ? &h_y : nullptr);

            // trajectory saving resets every generation: the first max_history_save
            // particles OF EACH GENERATION get their own entry, tagged with gen_idx
            if (save_history) {
                results.x_history.push_back(h_x);
                results.y_history.push_back(h_y);
                results.history_generation.push_back(gen_idx);
            }
        }

        print_progress(N_particles, gen_idx, true);

        // k estimate of this generation. The bank was produced with nu/k_used,
        // so bank/N alone measures k_true/k_used; multiplying by k_used undoes
        // the division and recovers k_true. (Using bank/N directly as the next k
        // makes k oscillate between k_true and ~1 instead of converging.)
        // It becomes the k the NEXT generation divides nu by.
        // (If the bank is empty the next iteration stops early before k_current
        // is ever used, so keep the old value instead of dividing by 0.)
        double k_est = k_current * static_cast<double>(fission_bank.size()) / N_particles;
        results.k_estimate.push_back(k_est);
        if (k_est > 0.0) k_current = k_est;
    }

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
