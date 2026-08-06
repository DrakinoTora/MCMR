#include "physics.hpp"
#include <cmath>
#include <random>

static std::random_device rd;
static std::mt19937 gen(rd());

double watt_spectrum_pdf(double E_MeV) {
    return 0.453 * std::exp(-1.036 * E_MeV) * std::sinh(std::sqrt(2.29 * E_MeV));
}

double watt_sample_single() {
    std::uniform_real_distribution<double> dist_E(0.0, 15.0);
    std::uniform_real_distribution<double> dist_Y(0.0, 0.4);
    while (true) {
        double E_rand = dist_E(gen);
        double Y_rand = dist_Y(gen);
        if (Y_rand <= watt_spectrum_pdf(E_rand)) {
            return E_rand * 1e6; // eV
        }
    }
}

double alpha_count(double A) {
    double term = (A - 1.0) / (A + 1.0);
    return term * term;
}

double E_scatter(double E, double alpha, double phi) {
    return E * ((1.0 + alpha) + (1.0 - alpha) * std::cos(phi)) / 2.0;
}

double Sigma_count(const std::vector<double>& E_data, const std::vector<double>& Sig_data, double E_target) {
    if (E_data.empty()) return 0.1;
    if (E_target <= E_data.front()) return Sig_data.front();
    if (E_target >= E_data.back()) return Sig_data.back();

    for (size_t i = 0; i < E_data.size() - 1; ++i) {
        if (E_target >= E_data[i] && E_target <= E_data[i+1]) {
            double t = (E_target - E_data[i]) / (E_data[i+1] - E_data[i]);
            return Sig_data[i] + t * (Sig_data[i+1] - Sig_data[i]);
        }
    }
    return Sig_data.front();
}