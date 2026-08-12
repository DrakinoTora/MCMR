#include "material.hpp"
#include <stdexcept>
#include <algorithm>

MaterialInfo get_material_info(const std::string& name) {
    std::string s = name;
    // Format to lowercase
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);

    if (s == "be" || s == "berilium") return {"Be", 425, 9.0};
    if (s == "c" || s == "grafit" || s == "carbon") return {"C", 600, 12.0};
    if (s == "fe" || s == "besi" || s == "iron") return {"Fe", 2631, 56.0};
    if (s == "pb" || s == "timbal" || s == "lead") return {"Pb", 8237, 208.0};

    throw std::invalid_argument("Material '" + name + "' unknown. Choose: Be, C, Fe, Pb");
}