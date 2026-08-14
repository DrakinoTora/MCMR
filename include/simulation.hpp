#pragma once
#include "particle.hpp"
#include "region.hpp"
#include "material.hpp"
#include "tally.hpp"
#include <map>
#include <string>
#include <vector>

class Simulation {
private:
    int N_particles;
    Grid grid;
    int max_history_save;

    // weight probability born neutron
    std::vector<double> flat_source_weights;

    std::map<int, std::vector<double>> E_data_total;
    std::map<int, std::vector<double>> Sig_data_total;
    std::map<int, std::vector<double>> E_data_scatter;
    std::map<int, std::vector<double>> Sig_data_scatter;

    Tally results;

public:
    Simulation(int N, double x_world, double y_world,
           const std::vector<double>& x_grid,
           const std::vector<double>& y_grid,
           const std::vector<std::vector<std::string>>& material_matrix,
           const std::vector<std::vector<double>>& sources,
           int max_save = 50,
           const std::string& bc_top   = "vacuum",
           const std::string& bc_bot   = "vacuum",
           const std::string& bc_left  = "vacuum",
           const std::string& bc_right = "vacuum");

    void set_cross_sections(
        const std::map<int, std::vector<double>>& E_tot,
        const std::map<int, std::vector<double>>& Sig_tot,
        const std::map<int, std::vector<double>>& E_scat,
        const std::map<int, std::vector<double>>& Sig_scat
    );

    void run();
    void export_xml(const std::string& filename = "mcmr_results.xml");

    Tally get_tally() const { return results; }
};
