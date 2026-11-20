#include "physics.hpp"
#include <algorithm>
#include <cmath>
#include <random>

double watt_spectrum_pdf(double E_MeV) {
    return 0.453 * std::exp(-1.036 * E_MeV) * std::sinh(std::sqrt(2.29 * E_MeV));
}

double watt_sample_single(std::mt19937& gen) {
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

void scatter_elastic(double A, double& E, double& mu_x, double& mu_y, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);

    const double mu_cm = 2.0 * dist_R(gen) - 1.0;           // isotropic in the CM frame
    const double alpha = alpha_count(A);
    E *= ((1.0 + alpha) + (1.0 - alpha) * mu_cm) / 2.0;     // energy after the collision

    // lab-frame deflection from the same CM angle
    const double denom = std::sqrt(A * A + 2.0 * A * mu_cm + 1.0);
    double cos_psi = (denom > 1e-12) ? (1.0 + A * mu_cm) / denom : 0.0;  // denom -> 0 only for A=1, mu_cm=-1
    cos_psi = std::clamp(cos_psi, -1.0, 1.0);
    double sin_psi = std::sqrt(std::max(0.0, 1.0 - cos_psi * cos_psi));
    if (dist_R(gen) < 0.5) sin_psi = -sin_psi;               // turn left or right of the incoming direction

    const double nx = mu_x * cos_psi - mu_y * sin_psi;
    const double ny = mu_x * sin_psi + mu_y * cos_psi;
    const double norm = std::sqrt(nx * nx + ny * ny);        // keep |mu| = 1 against rounding drift
    mu_x = nx / norm;
    mu_y = ny / norm;
}

int sample_fission_neutrons(double nu, std::mt19937& gen) {
    int base = static_cast<int>(std::floor(nu));
    double frac = nu - static_cast<double>(base);
    std::uniform_real_distribution<double> dist_R(0.0, 1.0);
    return (dist_R(gen) < frac) ? base + 1 : base;
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