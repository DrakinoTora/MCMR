#include "tally.hpp"
#include <algorithm>
#include <cmath>

void Tally::init_region_tally(int nx, int ny, double x_world, double y_world) {
    tally_nx = nx;
    tally_ny = ny;
    tally_dx = (nx > 0) ? x_world / nx : 0.0;
    tally_dy = (ny > 0) ? y_world / ny : 0.0;
    fission_tally.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), 0);
    absorp_tally.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), 0);
}

namespace {
size_t cell_index(double x, double y, int tally_nx, int tally_ny, double tally_dx, double tally_dy) {
    int ix_t = static_cast<int>(std::floor(x / tally_dx));
    int iy_t = static_cast<int>(std::floor(y / tally_dy));
    ix_t = std::clamp(ix_t, 0, tally_nx - 1);
    iy_t = std::clamp(iy_t, 0, tally_ny - 1);
    return static_cast<size_t>(iy_t) * tally_nx + ix_t;
}
}  // namespace

void Tally::add_fission_hit(double x, double y) {
    if (tally_nx <= 0 || tally_ny <= 0 || tally_dx <= 0.0 || tally_dy <= 0.0)
        return;
    fission_tally[cell_index(x, y, tally_nx, tally_ny, tally_dx, tally_dy)]++;
}

void Tally::add_absorp_hit(double x, double y) {
    if (tally_nx <= 0 || tally_ny <= 0 || tally_dx <= 0.0 || tally_dy <= 0.0)
        return;
    absorp_tally[cell_index(x, y, tally_nx, tally_ny, tally_dx, tally_dy)]++;
}