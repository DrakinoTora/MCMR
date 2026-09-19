#pragma once
#include "region.hpp"
#include "material.hpp"
#include "tally.hpp"
#include <deque>
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
//   - fission: the parent dies (tallied in fission_by_material) and
//     sample_fission_neutrons(nu[g]) children are pushed into fission_bank:
//     same position and same group as the parent, each with its own
//     isotropic direction. run() never starts the next source particle
//     while the bank is not empty -- it transports the banked neutrons
//     first (exactly like any other neutron; a fission there just appends
//     more members to the same bank). A supercritical setup therefore never
//     finishes, which is why run() prints the bank size next to the progress.
class SimulationMG {
private:
    // neutron born from fission, waiting in fission_bank to be transported
    struct BankedNeutron {
        double x, y;
        double mu_x, mu_y;
        int g;
    };

    int N_particles;
    Grid grid;
    int max_history_save;
    int n_groups = 0;

    std::vector<double> flat_source_weights;
    int n_grid_cells = 0;

    // per-material (by material symbol, e.g. "Fe") group data, each vector has n_groups elements
    std::map<std::string, std::vector<double>> sigma_t;
    std::map<std::string, std::vector<double>> sigma_s;
    std::map<std::string, std::vector<double>> sigma_f;
    std::map<std::string, std::vector<double>> nu;

    std::deque<BankedNeutron> fission_bank;  // FIFO queue

    Tally results;

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
           int max_save = 50,
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

    void run();
    void export_xml(const std::string& filename = "mcmr_results_mg.xml");

    Tally get_tally() const { return results; }
};
