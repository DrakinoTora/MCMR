#pragma once
#include "region.hpp"
#include "material.hpp"
#include "tally.hpp"
#include <map>
#include <random>
#include <string>
#include <vector>

// Multi-group (discretized energy) transport, parallel to Simulation (continuous
// energy). Shares Grid/Region/CircleRegion/Tally as-is -- only the energy
// treatment differs:
//   - particle born: group sampled uniform over [0, n_groups)
//   - sigma_t / sigma_s / sigma_f: direct per-group lookup, no interpolation
//   - scatter: isotropic direction, new group sampled uniform over [0, g]
//     (down-scatter only -- group 0 = fastest, higher index = lower energy)
//   - fission: the parent dies (tallied into the region tally, at its
//     position) and
//     sample_fission_neutrons(nu[g]) children are pushed into fission_bank:
//     same position and same group as the parent, each with its own
//     isotropic direction.
//
// GENERATIONS (power iteration / fission-source cycling, standard technique
// for k-eigenvalue problems): run(n_gen) processes n_gen generations, each
// with a FIXED population of N_particles neutrons.
//   - Generation 1: population sampled from the user-specified spatial/group
//     source distribution (exactly like the single-generation run used to).
//   - Generation g>1: population is N_particles draws WITH REPLACEMENT from
//     the fission_bank accumulated during generation g-1 (a bank member can
//     be picked more than once, or not at all) -- each draw keeps the bank
//     entry's birth position and group, but gets its own fresh isotropic
//     direction. This keeps the simulated population size constant across
//     generations regardless of whether the system is sub/super/critical.
//   - Within a generation, a fission does NOT get transported immediately --
//     it only appends children to fission_bank. All N_particles of the
//     current generation are transported first; the bank then becomes next
//     generation's source, and is cleared once population is drawn from it.
//   - If a generation's bank turns out empty (the whole population died out:
//     no fission anywhere), run() stops early -- there's nothing to sample
//     the next generation's population from.
//   - Trajectory saving (max_history_save) resets every generation: the
//     first max_history_save particles OF EACH GENERATION get a
//     <particle_history> entry, tagged with that generation's number.
class SimulationMG {
private:
    // a neutron ready to be transported: either one of this generation's
    // sampled source particles, or (before being resampled into the next
    // generation's population) a fission product waiting in fission_bank.
    struct BankedNeutron {
        double x, y;
        double mu_x, mu_y;
        int g;
    };

    int N_particles;
    Grid grid;
    int max_history_save;
    int n_groups = 0;

    // resolution of the region-based tally grid (see Tally::init_region_tally),
    // fixed at construction time, re-applied to `results` at the top of every run()
    int tally_res_x;
    int tally_res_y;

    std::vector<double> flat_source_weights;
    int n_grid_cells = 0;

    // per-material (by material symbol, e.g. "Fe") group data, each vector has n_groups elements
    std::map<std::string, std::vector<double>> sigma_t;
    std::map<std::string, std::vector<double>> sigma_s;
    std::map<std::string, std::vector<double>> sigma_f;
    std::map<std::string, std::vector<double>> nu;

    // fission products accumulated during the CURRENT generation, to be
    // resampled (with replacement) into next generation's population. No
    // longer a FIFO -- cleared and rebuilt fresh every generation.
    std::vector<BankedNeutron> fission_bank;

    Tally results;

    // k the CURRENT generation divides nu by (every fission samples
    // nu / k_current). run() sets it to its initial k for generation 1, then to
    // k_current * fission_bank.size() / N_particles at the end of each generation.
    double k_current = 1.0;

    // transport ONE neutron until it dies (leak / capture / fission).
    // h_x / h_y: trajectory is recorded into them, pass nullptr to not record.
    void transport_one(double x, double y, double mu_x, double mu_y, int g, int ix, int iy,
                       std::mt19937& gen,
                       std::vector<double>* h_x, std::vector<double>* h_y);

public:
    SimulationMG(int N, double x_world, double y_world,
           const std::vector<double>& x_grid,
           const std::vector<double>& y_grid,
           const std::vector<std::vector<std::string>>& material_matrix,
           const std::vector<std::vector<double>>& sources,
           const std::vector<double>& circle_cx = {},
           const std::vector<double>& circle_cy = {},
           const std::vector<double>& circle_r = {},
           const std::vector<std::string>& circle_material = {},
           const std::vector<double>& circle_source = {},
           int max_save = 10,
           int tally_res_x = 0,
           int tally_res_y = 0,
           const std::string& bc_top   = "vacuum",
           const std::string& bc_bot   = "vacuum",
           const std::string& bc_left  = "vacuum",
           const std::string& bc_right = "vacuum");

    void set_group_data(
        int n_groups_,
        const std::map<std::string, std::vector<double>>& sigma_t_,
        const std::map<std::string, std::vector<double>>& sigma_s_,
        const std::map<std::string, std::vector<double>>& sigma_f_,
        const std::map<std::string, std::vector<double>>& nu_
    );

    // n_gen: number of generations to cycle through (see the class comment
    // above for what a generation is). Default 1 = a single generation
    // sourced from the user-specified spatial/group distribution -- any
    // fission products are left sitting unresampled in fission_bank once
    // run() returns (their parent fission events are still reflected in the
    // region tally; the children themselves are NOT transported). Pass
    // n_gen>1 to cycle the fission
    // source across generations (power iteration).
    //
    // k: initial multiplication factor for generation 1 (default 1.0). Every
    // fission samples sample_fission_neutrons(nu / k_current); after each
    // generation k_current becomes k_current * fission_bank.size() /
    // N_particles for the next one (converges to the actual k). Each generation's k is recorded in Tally::k_used / k_estimate.
    void run(int n_gen = 1, double k = 1.0);
    void export_xml(const std::string& filename = "mcmr_results_mg.xml");

    Tally get_tally() const { return results; }
};
