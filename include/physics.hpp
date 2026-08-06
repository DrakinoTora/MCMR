#pragma once
#include <vector>

double watt_spectrum_pdf(double E_MeV);
double watt_sample_single();

double alpha_count(double A);
double E_scatter(double E, double alpha, double phi);
double Sigma_count(const std::vector<double>& E_data, const std::vector<double>& Sig_data, double E_target);