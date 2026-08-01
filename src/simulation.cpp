#include "simulation.hpp"
#include "geometry.hpp"
#include "physics.hpp"
#include <cmath>
#include <random>
#include <chrono>

SimulationResult neutron_sim(double a, double b, double c, int N_partikel, bool analog,
                            const std::vector<double>& E_Al_t, const std::vector<double>& Sig_Al_t,
                            const std::vector<double>& E_Al_s, const std::vector<double>& Sig_Al_s,
                            const std::vector<double>& E_Pb_t, const std::vector<double>& Sig_Pb_t,
                            const std::vector<double>& E_Pb_s, const std::vector<double>& Sig_Pb_s) 
{
    auto start_time = std::chrono::high_resolution_clock::now();
    SimulationResult res;

    res.E_start = watt_distribution_sampling(N_partikel);
    
    std::mt19937 gen(1337); // Seed tetap agar reproducible
    std::uniform_real_distribution<double> dist_x(b, a + b);
    std::uniform_real_distribution<double> dist_y(0.0, c);
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);

    for (int i = 0; i < N_partikel; ++i) {
        double w = 1.0;
        double E = res.E_start[i];
        double x = dist_x(gen);
        double y = dist_y(gen);

        double phi = dist_phi(gen);
        double mu_x = std::cos(phi);
        double mu_y = std::sin(phi);

        bool alive = true;
        int pos = 2; // Material A (Pb)

        while (alive) {
            double Sigma_t, Sigma_s;
            if (pos == 1 || pos == 3) {
                Sigma_t = Sigma_count(E_Al_t, Sig_Al_t, E);
                Sigma_s = Sigma_count(E_Al_s, Sig_Al_s, E);
            } else {
                Sigma_t = Sigma_count(E_Pb_t, Sig_Pb_t, E);
                Sigma_s = Sigma_count(E_Pb_s, Sig_Pb_s, E);
            }

            double R = dist_R(gen);
            double d = -std::log(R) / Sigma_t;

            x += d * mu_x;
            y += d * mu_y;

            auto [new_pos, A] = pos_particle(a, b, c, x, y);

            if (new_pos == 0) { // Kebocoran / Transmisi
                if (pos == 1 || pos == 3) res.E_B_escape.push_back(E);
                if (pos == 2) res.E_A_escape.push_back(E);
                res.transmisi += w;
                alive = false;
                break;
            }

            pos = new_pos;
            double P_scatter = Sigma_s / Sigma_t;

            if (analog) {
                if (dist_R(gen) <= P_scatter) {
                    double alpha = alpha_count(A);
                    phi = dist_phi(gen);
                    E = E_scatter(E, alpha, phi);
                    mu_x = std::cos(phi);
                    mu_y = std::sin(phi);
                } else {
                    if (pos == 2) res.abs_A += 1.0;
                    else res.abs_B += 1.0;
                    alive = false;
                }
            } else {
                w *= P_scatter;
                auto [status, new_w] = russian_roulette(w);
                w = new_w;

                if (status == "dead" || E < 0.025) {
                    alive = false;
                } else {
                    double alpha = alpha_count(A);
                    phi = dist_phi(gen);
                    E = E_scatter(E, alpha, phi);
                    mu_x = std::cos(phi);
                    mu_y = std::sin(phi);
                }

                double fraksi_abs = w * (1.0 - P_scatter);
                if (pos == 2) res.abs_A += fraksi_abs;
                else res.abs_B += fraksi_abs;
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    res.time_taken = diff.count();

    return res;
}