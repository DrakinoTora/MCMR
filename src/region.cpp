#include "region.hpp"
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cmath>

BoundaryType parse_boundary_type(const std::string& name) {
    std::string s = name;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "reflective" || s == "reflect" || s == "reflecting")
        return BoundaryType::Reflective;
    if (s == "vacuum" || s == "vacum" || s.empty())
        return BoundaryType::Vacuum;
    throw std::invalid_argument("Boundary condition '" + name + "' unknown. Choose: vacuum, reflective");
}

Grid::Grid(double x_world, double y_world,
           const std::vector<double>& x_grid,
           const std::vector<double>& y_grid,
           const std::vector<std::vector<std::string>>& material_matrix,
           const std::string& bc_left,
           const std::string& bc_right,
           const std::string& bc_bottom,
           const std::string& bc_top) {

    x_edges.push_back(0.0);
    x_edges.insert(x_edges.end(), x_grid.begin(), x_grid.end());
    x_edges.push_back(x_world);

    y_edges.push_back(0.0);
    y_edges.insert(y_edges.end(), y_grid.begin(), y_grid.end());
    y_edges.push_back(y_world);

    if (!std::is_sorted(x_edges.begin(), x_edges.end()))
        throw std::invalid_argument("x_grid must be in ascending order 0 to x_world");
    if (!std::is_sorted(y_edges.begin(), y_edges.end()))
        throw std::invalid_argument("y_grid must be in ascending order 0 to y_world");

    int nx = static_cast<int>(x_edges.size()) - 1;
    int ny = static_cast<int>(y_edges.size()) - 1;

    // material_matrix[row][col] convention:
    if (static_cast<int>(material_matrix.size()) != ny)
        throw std::invalid_argument(
            "material_matrix row count (" + std::to_string(material_matrix.size()) +
            ") must be same as len(y_grid)+1 = " + std::to_string(ny));

    regions.assign(nx, std::vector<Region>(ny));
    for (int row = 0; row < ny; ++row) {
        if (static_cast<int>(material_matrix[row].size()) != nx)
            throw std::invalid_argument(
                "material_matrix row " + std::to_string(row) +
                " col count must be same as len(x_grid)+1 = " + std::to_string(nx));

        int iy = ny - 1 - row; // top row (row=0) -> highest iy
        for (int col = 0; col < nx; ++col) {
            int ix = col; // left col (col=0) -> ix=0

            Region r;
            r.x1 = x_edges[ix];
            r.x2 = x_edges[ix + 1];
            r.y1 = y_edges[iy];
            r.y2 = y_edges[iy + 1];
            r.material = get_material_info(material_matrix[row][col]);
            regions[ix][iy] = r;
        }
    }

    bc_left_   = parse_boundary_type(bc_left);
    bc_right_  = parse_boundary_type(bc_right);
    bc_bottom_ = parse_boundary_type(bc_bottom);
    bc_top_    = parse_boundary_type(bc_top);
}

bool Grid::find_index(double x, double y, int& ix, int& iy) const {
    if (x < x_edges.front() || x > x_edges.back() ||
        y < y_edges.front() || y > y_edges.back()) {
        return false;
    }

    ix = static_cast<int>(std::upper_bound(x_edges.begin(), x_edges.end(), x) - x_edges.begin()) - 1;
    iy = static_cast<int>(std::upper_bound(y_edges.begin(), y_edges.end(), y) - y_edges.begin()) - 1;

    ix = std::clamp(ix, 0, static_cast<int>(regions.size()) - 1);
    iy = std::clamp(iy, 0, static_cast<int>(regions[0].size()) - 1);
    return true;
}

double Grid::distance_to_boundary(double x, double y, double mu_x, double mu_y,
                                   int ix, int iy, Side& hit_side) const {
    const Region& r = regions[ix][iy];
    double d_min = std::numeric_limits<double>::infinity();
    const double eps = 1e-12;
    hit_side = Side::None;

    const int nx = static_cast<int>(regions.size());
    const int ny = static_cast<int>(regions[0].size());

    if (mu_x > eps) {
        double d = (r.x2 - x) / mu_x;
        if (d < d_min) { d_min = d; hit_side = (ix == nx - 1) ? Side::Right : Side::None; }
    } else if (mu_x < -eps) {
        double d = (r.x1 - x) / mu_x;
        if (d < d_min) { d_min = d; hit_side = (ix == 0) ? Side::Left : Side::None; }
    }

    if (mu_y > eps) {
        double d = (r.y2 - y) / mu_y;
        if (d < d_min) { d_min = d; hit_side = (iy == ny - 1) ? Side::Top : Side::None; }
    } else if (mu_y < -eps) {
        double d = (r.y1 - y) / mu_y;
        if (d < d_min) { d_min = d; hit_side = (iy == 0) ? Side::Bottom : Side::None; }
    }

    // true circle boundary check: solve |(x,y) + t*(mu_x,mu_y) - (cx,cy)|^2 = r^2
    // for t (quadratic line-circle intersection), NOT a nearest-point distance.
    // (mu_x, mu_y) is a unit vector, so the quadratic's "a" coefficient is exactly 1.
    const double circ_eps = 1e-9;
    for (const auto& c : circles) {
        double dx = x - c.cx;
        double dy = y - c.cy;
        double b = 2.0 * (dx * mu_x + dy * mu_y);
        double cterm = dx * dx + dy * dy - c.r * c.r;
        double disc = b * b - 4.0 * cterm;
        if (disc < 0.0) continue;  // ray never touches this circle

        double sq = std::sqrt(disc);
        double t1 = (-b - sq) / 2.0;  // nearer intersection along the ray
        double t2 = (-b + sq) / 2.0;  // farther intersection along the ray

        // a circle boundary crossing is a material change, never a world edge
        if (t1 > circ_eps && t1 < d_min) { d_min = t1; hit_side = Side::None; }
        if (t2 > circ_eps && t2 < d_min) { d_min = t2; hit_side = Side::None; }
    }

    if (d_min <= 0.0) {
        hit_side = Side::None;
        return 1e-6;
    }
    return d_min;
}

void Grid::add_circle(double cx, double cy, double r, const std::string& material_name) {
    if (r <= 0.0)
        throw std::invalid_argument("circle radius must be positive");

    const double eps = 1e-9;
    if (cx - r < -eps || cx + r > world_max_x() + eps ||
        cy - r < -eps || cy + r > world_max_y() + eps) {
        throw std::invalid_argument(
            "circle centered at (" + std::to_string(cx) + ", " + std::to_string(cy) +
            ") with radius " + std::to_string(r) + " goes outside world bounds [0, " +
            std::to_string(world_max_x()) + "] x [0, " + std::to_string(world_max_y()) + "]");
    }

    CircleRegion c;
    c.cx = cx;
    c.cy = cy;
    c.r = r;
    c.material = get_material_info(material_name);
    circles.push_back(c);
}

const MaterialInfo& Grid::material_at(double x, double y, int ix, int iy) const {
    // check circles from last-registered to first -- later circles override earlier ones
    // wherever they overlap (painter's algorithm), same convention as region ops in Python
    for (int i = static_cast<int>(circles.size()) - 1; i >= 0; --i) {
        double dx = x - circles[i].cx;
        double dy = y - circles[i].cy;
        if (dx * dx + dy * dy <= circles[i].r * circles[i].r) {
            return circles[i].material;
        }
    }
    return regions[ix][iy].material;
}

BoundaryType Grid::bc_for_side(Side s) const {
    switch (s) {
        case Side::Left:   return bc_left_;
        case Side::Right:  return bc_right_;
        case Side::Bottom: return bc_bottom_;
        case Side::Top:    return bc_top_;
        default:           return BoundaryType::Vacuum;
    }
}