#include "region.hpp"
#include <algorithm>
#include <stdexcept>
#include <limits>

Grid::Grid(double x_world, double y_world,
           const std::vector<double>& x_grid,
           const std::vector<double>& y_grid,
           const std::vector<std::vector<std::string>>& material_matrix) {

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

    if (static_cast<int>(material_matrix.size()) != nx)
        throw std::invalid_argument(
            "material row (" + std::to_string(material_matrix.size()) +
            ") must be same as len(x_grid)+1 = " + std::to_string(nx));

    regions.resize(nx);
    for (int i = 0; i < nx; ++i) {
        if (static_cast<int>(material_matrix[i].size()) != ny)
            throw std::invalid_argument(
                "material col " + std::to_string(i) +
                " must be same as len(y_grid)+1 = " + std::to_string(ny));

        regions[i].resize(ny);
        for (int j = 0; j < ny; ++j) {
            Region r;
            r.x1 = x_edges[i];
            r.x2 = x_edges[i + 1];
            r.y1 = y_edges[j];
            r.y2 = y_edges[j + 1];
            r.material = get_material_info(material_matrix[i][j]);
            regions[i][j] = r;
        }
    }
}

bool Grid::find_index(double x, double y, int& ix, int& iy) const {
    if (x < x_edges.front() || x > x_edges.back() ||
        y < y_edges.front() || y > y_edges.back()) {
        return false; // outside world
    }

    ix = static_cast<int>(std::upper_bound(x_edges.begin(), x_edges.end(), x) - x_edges.begin()) - 1;
    iy = static_cast<int>(std::upper_bound(y_edges.begin(), y_edges.end(), y) - y_edges.begin()) - 1;

    // x == x_world still in cell
    ix = std::clamp(ix, 0, static_cast<int>(regions.size()) - 1);
    iy = std::clamp(iy, 0, static_cast<int>(regions[0].size()) - 1);
    return true;
}

double Grid::distance_to_boundary(double x, double y, double mu_x, double mu_y, int ix, int iy) const {
    const Region& r = regions[ix][iy];
    double d_min = std::numeric_limits<double>::infinity();
    const double eps = 1e-12;

    if (mu_x > eps)       d_min = std::min(d_min, (r.x2 - x) / mu_x);
    else if (mu_x < -eps) d_min = std::min(d_min, (r.x1 - x) / mu_x);

    if (mu_y > eps)       d_min = std::min(d_min, (r.y2 - y) / mu_y);
    else if (mu_y < -eps) d_min = std::min(d_min, (r.y1 - y) / mu_y);

    return (d_min <= 0.0) ? 1e-6 : d_min;
}
