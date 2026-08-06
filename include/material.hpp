#pragma once
#include <string>

struct MaterialInfo {
    std::string symbol;
    int mat_code;
    double A;
};

MaterialInfo get_material_info(const std::string& name);