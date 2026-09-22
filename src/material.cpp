#include "material.hpp"
#include <stdexcept>
#include <algorithm>

MaterialInfo get_material_info(const std::string& name) {
    std::string s = name;
    // Format to lowercase
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);

    // Atomic number density N = rho * N_A / M [atoms/cm^3], from standard
    // density (g/cm^3) and molar mass (g/mol). N_A = 6.02214076e23 /mol.
    //   Be: rho=1.848,  M=9.012182  -> N=1.234875e23
    //   C  (graphite): rho=2.260, M=12.0107   -> N=1.133159e23
    //   Fe: rho=7.874,  M=55.845   -> N=8.491062e22
    //   Pb: rho=11.34,  M=207.2    -> N=3.295901e22
    if (s == "be" || s == "berilium") return {"Be", 425, 9.0, 1.234875e23};
    if (s == "c" || s == "grafit" || s == "carbon") return {"C", 600, 12.0, 1.133159e23};
    if (s == "fe" || s == "besi" || s == "iron") return {"Fe", 2631, 56.0, 8.491062e22};
    if (s == "pb" || s == "timbal" || s == "lead") return {"Pb", 8237, 208.0, 3.295901e22};

    throw std::invalid_argument("Material '" + name + "' unknown. Choose: Be, C, Fe, Pb");
}