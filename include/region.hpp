#pragma once
#include "material.hpp"
#include <vector>
#include <string>

// region in grid: boundary + material
struct Region {
    double x1, y1, x2, y2;
    MaterialInfo material;
};

// true circular region overlaid on top of the rectangular grid -- NOT tied to
// any single grid cell, boundary is a real circle (checked via line-circle
// intersection), not an approximation
struct CircleRegion {
    double cx, cy, r;
    MaterialInfo material;
};

enum class BoundaryType { Vacuum, Reflective };
enum class Side { None, Left, Right, Bottom, Top };

BoundaryType parse_boundary_type(const std::string& name);

// Grid 2D
class Grid {
private:
    std::vector<double> x_edges;
    std::vector<double> y_edges;
    std::vector<std::vector<Region>> regions;
    std::vector<CircleRegion> circles;  // painter's algorithm: later entries override earlier ones

    BoundaryType bc_left_, bc_right_, bc_bottom_, bc_top_;

public:
    Grid(double x_world, double y_world,
         const std::vector<double>& x_grid,
         const std::vector<double>& y_grid,
         const std::vector<std::vector<std::string>>& material_matrix,
         const std::string& bc_left  = "vacuum",
         const std::string& bc_right = "vacuum",
         const std::string& bc_bottom = "vacuum",
         const std::string& bc_top   = "vacuum");

    double world_max_x() const { return x_edges.back(); }
    double world_max_y() const { return y_edges.back(); }
    int nx() const { return static_cast<int>(regions.size()); }
    int ny() const { return regions.empty() ? 0 : static_cast<int>(regions[0].size()); }

    bool find_index(double x, double y, int& ix, int& iy) const;

    const Region& region_at(int ix, int iy) const { return regions[ix][iy]; }

    // register a true circular region; throws if it doesn't fit entirely inside the world
    void add_circle(double cx, double cy, double r, const std::string& material_name);
    int num_circles() const { return static_cast<int>(circles.size()); }
    const CircleRegion& circle_at_index(int i) const { return circles[i]; }

    // actual material at a point, accounting for circles drawn on top of the grid cell
    const MaterialInfo& material_at(double x, double y, int ix, int iy) const;

    // distance to boundary; hit_side is set when the boundary hit is a world edge.
    // Also checks every circle for a line-circle intersection along (mu_x, mu_y) --
    // whichever event (rectangle edge or circle boundary) is closer wins.
    double distance_to_boundary(double x, double y, double mu_x, double mu_y,
                                 int ix, int iy, Side& hit_side) const;

    BoundaryType bc_for_side(Side s) const;
};