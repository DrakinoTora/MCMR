#pragma once
#include <random>
#include <vector>

double watt_spectrum_pdf(double E_MeV);
// draws from the CALLER's generator (the run's seeded RNG), never a private one
double watt_sample_single(std::mt19937& gen);

double alpha_count(double A);
double E_scatter(double E, double alpha, double phi);
double Sigma_count(const std::vector<double>& E_data, const std::vector<double>& Sig_data, double E_target);

// Stochastic rounding of nu (mean neutrons per fission) so the expected value
int sample_fission_neutrons(double nu, std::mt19937& gen);