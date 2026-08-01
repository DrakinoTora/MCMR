#pragma once
#include "particle.hpp"

// Fungsi utama simulasi Monte Carlo
SimulationResult neutron_sim(double a, double b, double c, int N_partikel, bool analog,
                            const std::vector<double>& E_Al_t, const std::vector<double>& Sig_Al_t,
                            const std::vector<double>& E_Al_s, const std::vector<double>& Sig_Al_s,
                            const std::vector<double>& E_Pb_t, const std::vector<double>& Sig_Pb_t,
                            const std::vector<double>& E_Pb_s, const std::vector<double>& Sig_Pb_s);