#include "tally.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

void Tally::init_region_tally(int nx, int ny, double x_world, double y_world) {
    tally_nx = nx;
    tally_ny = ny;
    tally_dx = (nx > 0) ? x_world / nx : 0.0;
    tally_dy = (ny > 0) ? y_world / ny : 0.0;
    fission_tally.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), 0);
    absorp_tally.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), 0);
    flux_tally.assign(static_cast<size_t>(nx) * static_cast<size_t>(ny), 0.0);
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

void Tally::add_flux_track(double x0, double y0, double x1, double y1) {
    if (tally_nx <= 0 || tally_ny <= 0 || tally_dx <= 0.0 || tally_dy <= 0.0)
        return;

    const double x_world = tally_dx * tally_nx;
    const double y_world = tally_dy * tally_ny;

    // clamp both endpoints onto the grid first -- a neutron sitting exactly on
    // (or, from floating point, a hair outside) the world edge must not
    // desync the traversal below
    x0 = std::clamp(x0, 0.0, x_world); y0 = std::clamp(y0, 0.0, y_world);
    x1 = std::clamp(x1, 0.0, x_world); y1 = std::clamp(y1, 0.0, y_world);

    const double dx = x1 - x0;
    const double dy = y1 - y0;
    const double length = std::sqrt(dx * dx + dy * dy);
    if (length < 1e-12) return;  // no movement -- nothing to score

    // Amanatides & Woo voxel traversal: walk the segment cell by cell in
    // parametric t in [0,1], scoring each cell's exact overlap length.
    int ix = std::clamp(static_cast<int>(std::floor(x0 / tally_dx)), 0, tally_nx - 1);
    int iy = std::clamp(static_cast<int>(std::floor(y0 / tally_dy)), 0, tally_ny - 1);

    const int step_x = (dx > 0.0) ? 1 : (dx < 0.0 ? -1 : 0);
    const int step_y = (dy > 0.0) ? 1 : (dy < 0.0 ? -1 : 0);

    // t_max_* = t at which the segment next crosses a vertical (resp.
    // horizontal) grid line; t_delta_* = t needed to cross one full cell.
    // An axis with zero movement never crosses a line -> stays at +inf.
    constexpr double kInf = std::numeric_limits<double>::infinity();
    double t_max_x = kInf, t_delta_x = kInf;
    if (step_x != 0) {
        double next_x = (ix + (step_x > 0 ? 1 : 0)) * tally_dx;
        t_max_x = (next_x - x0) / dx;
        t_delta_x = tally_dx / std::fabs(dx);
    }
    double t_max_y = kInf, t_delta_y = kInf;
    if (step_y != 0) {
        double next_y = (iy + (step_y > 0 ? 1 : 0)) * tally_dy;
        t_max_y = (next_y - y0) / dy;
        t_delta_y = tally_dy / std::fabs(dy);
    }

    double t = 0.0;
    while (t < 1.0 - 1e-12) {
        double t_next = std::min({t_max_x, t_max_y, 1.0});

        double seg_len = (t_next - t) * length;
        if (seg_len > 0.0)
            flux_tally[static_cast<size_t>(iy) * tally_nx + ix] += seg_len;

        t = t_next;
        // advance whichever axis (or both, on an exact corner hit) hit its line
        bool cross_x = (t_max_x <= t_next + 1e-12);
        bool cross_y = (t_max_y <= t_next + 1e-12);
        if (cross_x) { ix += step_x; t_max_x += t_delta_x; }
        if (cross_y) { iy += step_y; t_max_y += t_delta_y; }

        // left the grid -- a straight segment can't re-enter
        if (ix < 0 || ix >= tally_nx || iy < 0 || iy >= tally_ny) break;
    }
}
