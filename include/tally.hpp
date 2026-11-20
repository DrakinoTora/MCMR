#pragma once
#include <cstdint>
#include <vector>
#include <string>

struct Tally {
    // World leakage (vacuum boundaries only), split by which side of the world
    // the neutron escaped through. Replaces the old single `transmission`
    // counter -- the total leak is leak_left + leak_right + leak_bottom + leak_top
    // (= E_leak.size() in continuous mode, G_leak.size() in group mode).
    int leak_left   = 0;  // escaped through x = 0
    int leak_right  = 0;  // escaped through x = x_world
    int leak_bottom = 0;  // escaped through y = 0
    int leak_top    = 0;  // escaped through y = y_world
    double time_taken = 0.0;
    std::uint64_t seed = 0; // seed the run actually used (given, or drawn at random) -> pass it back as seed= to repeat the run exactly

    std::vector<double> E_born;  // continuous mode only -- birth energy (eV) per source particle
    std::vector<double> E_leak;  // continuous mode only -- leak energy (eV), one entry per leaked neutron (all 4 sides)

    std::vector<int> G_born;     // group mode only -- birth group index per source particle
    std::vector<int> G_leak;     // group mode only -- leak group index, one entry per leaked neutron (all 4 sides)

    std::vector<std::vector<double>> x_history;
    std::vector<std::vector<double>> y_history;
    // parallel to x_history/y_history: which generation (1-based) that
    // particle_history entry belongs to. Group mode only -- always empty
    // for continuous mode (single generation, no concept of "gen").
    std::vector<int> history_generation;

    // Group mode only (empty for continuous): one entry per generation that
    // actually ran.
    //   k_used[i]     = the k generation i+1 divided nu by (gen 1: run()'s
    //                   initial k; gen g>1: k_estimate of gen g-1)
    //   k_estimate[i] = k_used[i] * fission_bank size at the end of
    //                   generation i+1 / N_particles  (the bank was made with
    //                   nu/k_used, so this recovers the actual k)
    std::vector<double> k_used;
    std::vector<double> k_estimate;

    // ------------------------------------------------------------------ //
    // Region-based tally grids: a resolution grid laid independently over the
    // WHOLE WORLD (not tied to material_matrix's grid at all). Fission and
    // absorption events are tracked in SEPARATE grids so each can be
    // visualized (heatmap) independently -- fission_tally only counts
    // fission events, absorp_tally only counts (non-fission) absorption
    // events. Continuous mode (Simulation, no fission model) only ever
    // populates absorp_tally and flux_tally; fission_tally stays all-zero
    // for it.
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
    // Track-length flux estimator, PER ENERGY GROUP, on the same cell grid as
    // fission_tally/absorp_tally (same tally_nx/tally_ny/tally_dx/tally_dy).
    // Every straight-line segment a neutron travels (birth-to-collision,
    // collision-to-collision, collision-to-boundary, reflection-to-collision,
    // ...) adds its length to EVERY cell that segment passes through, weighted
    // by exactly how much of the segment lies in that cell, and it is scored
    // into the grid of the group the neutron is IN while travelling that
    // segment (a neutron only changes group at a scatter, i.e. between
    // segments). E.g. a group-2 neutron crossing a cell for 2 cm adds 2.0 to
    // that cell in group 2's grid.
    //
    // Flat layout: index = g * (tally_nx*tally_ny) + iy_t * tally_nx + ix_t
    // (g = 0 is the fastest group, same numbering as G_born/G_leak), size =
    // flux_groups * tally_nx * tally_ny. Group mode: flux_groups = n_groups.
    // Continuous mode has no groups: flux_groups = 1, a single grid holding
    // the total flux. Summing over g gives the total (all-group) flux.
    std::vector<double> flux_tally;
    int flux_groups = 1;  // number of group grids stored in flux_tally (and the two rate tallies below)

    // Reaction-RATE tallies (track-length estimator of  \int Sigma_x * phi dV):
    // same layout, groups and cell grid as flux_tally, but each piece of a
    // segment scores  (length inside the cell) * Sigma_x  instead of just the
    // length, where Sigma_x is the macroscopic cross section [1/cm] of the
    // material the segment is travelling through, for the neutron's current
    // energy/group:
    //   fission_rate_tally : Sigma_f (group mode; stays all-zero in continuous
    //                        mode, which has no fission model)
    //   absorp_rate_tally  : Sigma_a = pure capture, i.e. what the engine itself
    //                        treats as non-fission absorption -- group mode:
    //                        Sigma_t - Sigma_s - Sigma_f; continuous mode:
    //                        Sigma_t(E) - Sigma_s(E)  (clamped at >= 0)
    // A segment never spans two materials (the engine stops a neutron at every
    // material cell / circle boundary), so Sigma is constant along a segment; when
    // a neutron crosses into another material, the NEXT segment simply uses
    // that material's Sigma. Unlike fission_tally/absorp_tally (integer counts
    // of events), these score wherever a neutron merely passes through, so they
    // are far less noisy for the same N.
    std::vector<double> fission_rate_tally;
    std::vector<double> absorp_rate_tally;
    int tally_nx = 0;
    int tally_ny = 0;
    double tally_dx = 0.0;
    double tally_dy = 0.0;

    // (Re)initialize all tally grids to nx x ny cells covering [0,x_world] x
    // [0,y_world], all zeroed. flux_tally gets n_groups grids (default 1). Called once per run() (after the rest of the
    // Tally is reset) so every fresh run starts from an all-zero grid.
    void init_region_tally(int nx, int ny, double x_world, double y_world, int n_groups = 1);

    // Increment the fission_tally / absorp_tally cell containing (x, y) by 1.
    // A point that falls (very slightly, from floating point) outside
    // [0,x_world] x [0,y_world] is clamped to the nearest edge cell instead
    // of being silently dropped.
    void add_fission_hit(double x, double y);
    void add_absorp_hit(double x, double y);

    // Track-length score of the straight segment (x0,y0)->(x1,y1) in `group`.
    // Per cell crossed (with L = length of the segment inside that cell):
    //   flux_tally          += L
    //   fission_rate_tally  += L * sigma_f
    //   absorp_rate_tally   += L * sigma_a
    // (one grid traversal scores all three). Original description:
    // Track-length flux score: add the length of the straight segment from
    // (x0,y0) to (x1,y1) to every cell of flux_tally's grid for `group` it crosses, split
    // proportionally to how much of the segment lies in each cell. A
    // zero-length segment is a no-op. Endpoints outside
    // [0,x_world] x [0,y_world] are clamped onto the grid first.
    // A group outside [0, flux_groups) is ignored.
    void add_flux_track(double x0, double y0, double x1, double y1, int group = 0,
                        double sigma_f = 0.0, double sigma_a = 0.0);
};