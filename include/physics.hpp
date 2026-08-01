#pragma once
#include <vector>
#include <string>

// Pembangkitan energi Watt Spectrum
double watt_spectrum_pdf(double E_MeV);
std::vector<double> watt_distribution_sampling(int neutron_counts);

// Perhitungan fisika
double alpha_count(double A);
double E_scatter(double E, double alpha, double phi);
std::pair<std::string, double> russian_roulette(double weight, double w_cutoff = 0.001, double p_kill = 0.5);

// Interpolasi Cross-Section (pengganti Sigma_count)
double Sigma_count(const std::vector<double>& E_data, const std::vector<double>& Sig_data, double E_target);