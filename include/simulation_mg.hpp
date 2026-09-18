#pragma once
#include "region.hpp"
#include "material.hpp"
#include "tally.hpp"
#include <map>
#include <string>
#include <vector>

// Multi-group (discretized energy) transport, parallel to Simulation (continuous
// energy). Shares Grid/Region/CircleRegion/Tally as-is -- only the energy
// treatment differs:
//   - particle born: group sampled uniform over [0, n_groups)
//   - sigma_t / sigma_s / sigma_f: direct per-group lookup, no interpolation
//   - scatter: isotropic direction, new group sampled uniform over [0, g]
//     (down-scatter only -- group 0 = fastest, higher index = lower energy)
//   - fission: sigma_f only subtracts from sigma_a for now (particle is
//     absorbed there, tallied separately in fission_by_material); no neutron
//     multiplication yet -- sample_fission_neutrons() in physics.hpp is ready
//     for when that's switched on.
class SimulationMG {
private:
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

    Tally results;

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
