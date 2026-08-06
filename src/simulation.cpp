#include "simulation.hpp"
#include "physics.hpp"
#include "exporter.hpp"
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>
#include <sstream>

#include <pybind11/pybind11.h>
namespace py = pybind11;

Simulation::Simulation(int N, const std::string& mat_name, BoxBoundary b, WorldBoundary w, int max_save)
    : N_particles(N), box(b), world(w), max_history_save(max_save) {
    mat_box = get_material_info(mat_name);
    mat_outside = get_material_info("Pb");
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
    std::uniform_real_distribution<double> dist_x(box.x1, box.x2);
    std::uniform_real_distribution<double> dist_y(box.y1, box.y2);
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);

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

        double x = dist_x(gen);
        double y = dist_y(gen);
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
            bool in_box = is_inside_box(x, y, box);
            MaterialInfo cur_mat = in_box ? mat_box : mat_outside;

            double Sigma_t = Sigma_count(E_data_total.at(cur_mat.mat_code), Sig_data_total.at(cur_mat.mat_code), E);
            double Sigma_s = Sigma_count(E_data_scatter.at(cur_mat.mat_code), Sig_data_scatter.at(cur_mat.mat_code), E);

            double R = dist_R(gen);
            double d_coll = -std::log(R) / Sigma_t;
            double d_surf = distance_to_surface(x, y, mu_x, mu_y, box, world);

            double d = std::min(d_coll, d_surf);
            x += d * mu_x;
            y += d * mu_y;

            if (is_outside_world(x, y, world)) {
                results.transmisi++;
                results.E_leak.push_back(E);
                alive = false;
                break;
            }

            if (d_coll >= d_surf) {
                x += 1e-5 * mu_x;
                y += 1e-5 * mu_y;
                continue;
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
                if (in_box) results.absorp_material++;
                else results.absorp_outside++;
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