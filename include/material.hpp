#pragma once
#include <string>

struct MaterialInfo {
    std::string symbol;
    // Key into the continuous-energy cross-section data (data/<mat_code>-N.html).
    // -1 means this element has no such data: it is still a valid material for
    // multi-group mode, but can't be used in continuous-energy mode.
    int mat_code;
    double A;
    // Atomic number density N [atoms/cm^3], used to convert the microscopic
    // cross sections loaded from data (barn) into macroscopic Sigma (cm^-1):
    // Sigma[cm^-1] = sigma[barn] * atom_density * 1e-24.
    double atom_density;
};

MaterialInfo get_material_info(const std::string& name);