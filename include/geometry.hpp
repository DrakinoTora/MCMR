#pragma once
#include "particle.hpp"

struct BoxBoundary {
    double x1, y1, x2, y2;
};

struct WorldBoundary {
    double max_x, max_y;
};

bool is_inside_box(double x, double y, const BoxBoundary& box);
bool is_outside_world(double x, double y, const WorldBoundary& world);

double distance_to_surface(double x, double y, double mu_x, double mu_y, 
                           const BoxBoundary& box, const WorldBoundary& world);