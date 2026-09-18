#pragma once

struct Particle {
    double x, y;
    double mu_x, mu_y;
    double E; // eV
    int g;    // energy group index
    bool alive;
};