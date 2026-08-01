#include "geometry.hpp"

std::pair<int, double> pos_particle(double a, double b, double c, double x, double y) {
    if (y >= 0.0 && y <= c) {
        if (x >= b && x <= (b + a)) {
            return {2, 208.0}; // Pb (Lead)
        } else if (x >= 0.0 && x <= b) {
            return {1, 56.0};  // Fe (Iron) Kiri
        } else if (x >= (b + a) && x <= (2 * b + a)) {
            return {3, 56.0};  // Fe (Iron) Kanan
        }
    }
    return {0, 0.0}; // Bocor / Luar Medium
}