#pragma once
#include <vector>

struct Tally {
    int absorp_material = 0;
    int absorp_outside = 0;
    int transmisi = 0;
    double time_taken = 0.0;

    std::vector<double> E_born;
    std::vector<double> E_leak;

    std::vector<std::vector<double>> x_history;
    std::vector<std::vector<double>> y_history;
};