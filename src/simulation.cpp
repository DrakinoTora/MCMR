#include "simulation.hpp"
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

Simulation::Simulation(int N, double x_world, double y_world,
                        const std::vector<double>& x_grid,
                        const std::vector<double>& y_grid,
                        const std::vector<std::vector<std::string>>& material_matrix,
                        const std::vector<std::vector<double>>& sources,
                        int max_save,
                        const std::string& bc_top,
                        const std::string& bc_bot,
                        const std::string& bc_left,
                        const std::string& bc_right)
    : N_particles(N),
      grid(x_world, y_world, x_grid, y_grid, material_matrix,
           bc_left, bc_right, bc_bot, bc_top),
      max_history_save(max_save) {

    int nx = grid.nx();
    int ny = grid.ny();

    if (static_cast<int>(sources.size()) != ny)
        throw std::invalid_argument(
            "sources row count (" + std::to_string(sources.size()) +
            ") must be same as len(y_grid)+1 = " + std::to_string(ny));

    flat_source_weights.resize(nx * ny);
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
    if (total <= 0.0)
        throw std::invalid_argument("sum all sources must be greater than 0");
}

void Simulation::set_cross_sections(
    const std::map<int, std::vector<double>>& E_tot,
    const std::map<int, std::vector<double>>& Sig_tot,
    const std::map<int, std::vector<double>>& E_scat,
    const std::map<int, std::vector<double>>& Sig_scat
) {
    E_data_total = E_tot;
    Sig_data_total = Sig_tot;
    E_data_scatter = E_scat;
    Sig_data_scatter = Sig_scat;
}

void Simulation::run() {
    py::print("====================================================");
    py::print("              MCMR Simulation Engine               ");
    py::print("====================================================");
    py::print(" Total Particle   :", N_particles);
    py::print("====================================================");
    py::print("running...");

    py::module_::import("sys").attr("stdout").attr("flush")();

    auto start_time = std::chrono::high_resolution_clock::now();

    std::mt19937 gen(1337);
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);

    // born region
    std::discrete_distribution<int> region_picker(flat_source_weights.begin(), flat_source_weights.end());
    const int ny = grid.ny();

    int step_update = N_particles / 100;
    if (step_update < 1) step_update = 1;

    for (int i = 0; i < N_particles; ++i) {

        if ((i + 1) % step_update == 0 || i == N_particles - 1) {
            int current = i + 1;
            int percent = (current * 100) / N_particles;

            std::stringstream ss;
            ss << "\rProgress: [" << current << "/" << N_particles << "] (" << percent << "%)";

            py::print(ss.str(), py::arg("end") = "");
            py::module_::import("sys").attr("stdout").attr("flush")();
        }

        int flat_idx = region_picker(gen);
        int ix = flat_idx / ny;
        int iy = flat_idx % ny;
        const Region& birth_region = grid.region_at(ix, iy);

        std::uniform_real_distribution<double> dist_bx(birth_region.x1, birth_region.x2);
        std::uniform_real_distribution<double> dist_by(birth_region.y1, birth_region.y2);
        double x = dist_bx(gen);
        double y = dist_by(gen);

        double phi = dist_phi(gen);
        double mu_x = std::cos(phi);
        double mu_y = std::sin(phi);
        double E = watt_sample_single();

        results.E_born.push_back(E);

        std::vector<double> h_x, h_y;
        bool save_history = (i < max_history_save);
        if (save_history) {
            h_x.push_back(x);
            h_y.push_back(y);
        }

        bool alive = true;
        while (alive) {
            const Region& cur_region = grid.region_at(ix, iy);
            const MaterialInfo& cur_mat = cur_region.material;

            double Sigma_t = Sigma_count(E_data_total.at(cur_mat.mat_code), Sig_data_total.at(cur_mat.mat_code), E);
            double Sigma_s = Sigma_count(E_data_scatter.at(cur_mat.mat_code), Sig_data_scatter.at(cur_mat.mat_code), E);

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
                        results.E_leak.push_back(E);
                        alive = false;
                        break;
                    } else {
                        if (hit_side == Side::Left || hit_side == Side::Right) {
                            mu_x = -mu_x;
                        } else {
                            mu_y = -mu_y;
                        }
                        if (save_history) {
                            h_x.push_back(x);
                            h_y.push_back(y);
                        }
                        continue;
                    }
                }

                x += 1e-6 * mu_x;
                y += 1e-6 * mu_y;
                grid.find_index(x, y, ix, iy);
                continue;
            }

            if (save_history) {
                h_x.push_back(x);
                h_y.push_back(y);
            }

            if (save_history) {
                h_x.push_back(x);
                h_y.push_back(y);
            }

            double P_scatter = Sigma_s / Sigma_t;
            double P = dist_R(gen);

            if (P < P_scatter) {
                phi = dist_phi(gen);
                mu_x = std::cos(phi);
                mu_y = std::sin(phi);
                double alpha = alpha_count(cur_mat.A);
                E = E_scatter(E, alpha, phi);
            } else {
                alive = false;
                results.absorp_by_material[cur_mat.symbol]++;
            }
        }

        if (save_history) {
            results.x_history.push_back(h_x);
            results.y_history.push_back(h_y);
        }
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

    export_xml();
}

void Simulation::export_xml(const std::string& filename) {
    export_to_xml(results, filename);
}
