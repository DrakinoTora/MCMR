#pragma once
#include <vector>
#include <map>
#include <string>

struct Tally {
    // neutron absorb per material
    std::map<std::string, int> absorp_by_material;
    int transmisi = 0; // world leak
    double time_taken = 0.0;

    std::vector<double> E_born;
    std::vector<double> E_leak;

    std::vector<std::vector<double>> x_history;
    std::vector<std::vector<double>> y_history;
};