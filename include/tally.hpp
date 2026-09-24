#pragma once
#include <vector>
#include <map>
#include <string>

struct Tally {
    // neutron absorb per material
    std::map<std::string, int> absorp_by_material;
    std::map<std::string, int> fission_by_material;
    int transmission = 0; // world leak
    double time_taken = 0.0;

    std::vector<double> E_born;  // continuous mode only -- birth energy (eV) per source particle
    std::vector<double> E_leak;  // continuous mode only -- leak energy (eV), len == transmission

    std::vector<int> G_born;     // group mode only -- birth group index per source particle
    std::vector<int> G_leak;     // group mode only -- leak group index, len == transmission

    std::vector<std::vector<double>> x_history;
    std::vector<std::vector<double>> y_history;
};