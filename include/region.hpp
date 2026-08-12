#pragma once
#include "material.hpp"
#include <vector>
#include <string>

// region in grid: boundary + material
struct Region {
    double x1, y1, x2, y2;
    MaterialInfo material;
};

// Grid 2D
class Grid {
private:
    std::vector<double> x_edges; // x grid include 0 and world edge
    std::vector<double> y_edges; // y grid include 0 and world edge
    std::vector<std::vector<Region>> regions; // regions[ix][iy]

public:
    Grid(double x_world, double y_world,
         const std::vector<double>& x_grid,
         const std::vector<double>& y_grid,
         const std::vector<std::vector<std::string>>& material_matrix);

    double world_max_x() const { return x_edges.back(); }
    double world_max_y() const { return y_edges.back(); }

    // index region
    bool find_index(double x, double y, int& ix, int& iy) const;

    const Region& region_at(int ix, int iy) const { return regions[ix][iy]; }

    // d surface
    double distance_to_boundary(double x, double y, double mu_x, double mu_y, int ix, int iy) const;
};
