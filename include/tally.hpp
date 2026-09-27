#pragma once
#include <vector>
#include <string>

struct Tally {
    int transmission = 0; // world leak
    double time_taken = 0.0;

    std::vector<double> E_born;  // continuous mode only -- birth energy (eV) per source particle
    std::vector<double> E_leak;  // continuous mode only -- leak energy (eV), len == transmission

    std::vector<int> G_born;     // group mode only -- birth group index per source particle
    std::vector<int> G_leak;     // group mode only -- leak group index, len == transmission

    std::vector<std::vector<double>> x_history;
    std::vector<std::vector<double>> y_history;
    // parallel to x_history/y_history: which generation (1-based) that
    // particle_history entry belongs to. Group mode only -- always empty
    // for continuous mode (single generation, no concept of "gen").
    std::vector<int> history_generation;

    // ------------------------------------------------------------------ //
    // Region-based tally grids: a resolution grid laid independently over the
    // WHOLE WORLD (not tied to material_matrix's grid at all). Fission and
    // absorption events are tracked in SEPARATE grids so each can be
    // visualized (heatmap) independently -- fission_tally only counts
    // fission events, absorp_tally only counts (non-fission) absorption
    // events. Continuous mode (Simulation, no fission model) only ever
    // populates absorp_tally; fission_tally stays all-zero for it.
    // tally_nx/tally_ny = number of columns/rows; tally_dx/tally_dy = each
    // cell's width/height in world units (= x_world/tally_nx, y_world/tally_ny).
    //
    // Flat, row-major: index = iy_t * tally_nx + ix_t, where ix_t=0 is the
    // LEFTMOST column (x=0) and iy_t=0 is the BOTTOM row (y=0) -- this
    // matches the engine's own internal ix/iy convention (NOT the
    // row=0=topmost convention used by material_matrix/sources in the
    // Python API).
    // ------------------------------------------------------------------ //
    std::vector<int> fission_tally;
    std::vector<int> absorp_tally;
    int tally_nx = 0;
    int tally_ny = 0;
    double tally_dx = 0.0;
    double tally_dy = 0.0;

    // (Re)initialize both tally grids to nx x ny cells covering [0,x_world] x
    // [0,y_world], all zeroed. Called once per run() (after the rest of the
    // Tally is reset) so every fresh run starts from an all-zero grid.
    void init_region_tally(int nx, int ny, double x_world, double y_world);

    // Increment the fission_tally / absorp_tally cell containing (x, y) by 1.
    // A point that falls (very slightly, from floating point) outside
    // [0,x_world] x [0,y_world] is clamped to the nearest edge cell instead
    // of being silently dropped.
    void add_fission_hit(double x, double y);
    void add_absorp_hit(double x, double y);
};