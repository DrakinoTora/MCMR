#pragma once
#include <vector>

struct Particle {
    double x, y;
    double mu_x, mu_y;
    double E;          // Energi dalam eV
    double weight;     // Bobot neutron
    bool alive;
    int pos;           // 1: Fe Kiri, 2: Pb, 3: Fe Kanan, 0: Bocor/Lolos
};

struct SimulationResult {
    double abs_A = 0.0;
    double abs_B = 0.0;
    double transmisi = 0.0;
    double time_taken = 0.0;
    std::vector<double> E_start;
    std::vector<double> E_A_escape;
    std::vector<double> E_B_escape;
};