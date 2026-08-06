#include "geometry.hpp"
#include <algorithm>
#include <cmath>

bool is_inside_box(double x, double y, const BoxBoundary& box) {
    return (x >= box.x1 && x <= box.x2 && y >= box.y1 && y <= box.y2);
}

bool is_outside_world(double x, double y, const WorldBoundary& world) {
    return (x < 0.0 || x > world.max_x || y < 0.0 || y > world.max_y);
}

double distance_to_surface(double x, double y, double mu_x, double mu_y, 
                           const BoxBoundary& box, const WorldBoundary& world) {
    double d_min = 1e9;

    if (mu_x > 0) {
        if (x < box.x1) d_min = std::min(d_min, (box.x1 - x) / mu_x);
        if (x < box.x2) d_min = std::min(d_min, (box.x2 - x) / mu_x);
        d_min = std::min(d_min, (world.max_x - x) / mu_x);
    } else if (mu_x < 0) {
        if (x > box.x2) d_min = std::min(d_min, (box.x2 - x) / mu_x);
        if (x > box.x1) d_min = std::min(d_min, (box.x1 - x) / mu_x);
        d_min = std::min(d_min, (0.0 - x) / mu_x);
    }

    if (mu_y > 0) {
        if (y < box.y1) d_min = std::min(d_min, (box.y1 - y) / mu_y);
        if (y < box.y2) d_min = std::min(d_min, (box.y2 - y) / mu_y);
        d_min = std::min(d_min, (world.max_y - y) / mu_y);
    } else if (mu_y < 0) {
        if (y > box.y2) d_min = std::min(d_min, (box.y2 - y) / mu_y);
        if (y > box.y1) d_min = std::min(d_min, (box.y1 - y) / mu_y);
        d_min = std::min(d_min, (0.0 - y) / mu_y);
    }

    return (d_min <= 0) ? 1e-5 : d_min;
}