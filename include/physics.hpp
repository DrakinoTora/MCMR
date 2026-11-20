#pragma once
#include <random>
#include <vector>

double watt_spectrum_pdf(double E_MeV);
// draws from the CALLER's generator (the run's seeded RNG), never a private one
double watt_sample_single(std::mt19937& gen);

double alpha_count(double A);
double E_scatter(double E, double alpha, double phi);

// Elastic scattering off a nucleus of mass number A (target at rest), isotropic
// in the centre-of-mass frame. ONE CM-frame draw decides BOTH outcomes, so
// energy and direction stay physically tied together:
//   * mu_cm = cos(theta_cm), uniform in [-1, 1]
//   * E'/E  = (A^2 + 2 A mu_cm + 1) / (A+1)^2  = ((1+alpha) + (1-alpha) mu_cm) / 2
//   * lab deflection angle psi from the INCOMING direction:
//       cos(psi) = (1 + A mu_cm) / sqrt(A^2 + 2 A mu_cm + 1)
// (mu_x, mu_y) is rotated by +-psi (side chosen at random: the engine tracks
// direction in the x-y plane only). Nothing depends on the absolute direction
// of travel, so a neutron heading +x and one heading -x lose energy alike.
void scatter_elastic(double A, double& E, double& mu_x, double& mu_y, std::mt19937& gen);
double Sigma_count(const std::vector<double>& E_data, const std::vector<double>& Sig_data, double E_target);

// Stochastic rounding of nu (mean neutrons per fission) so the expected value
int sample_fission_neutrons(double nu, std::mt19937& gen);