#pragma once
#include "particle.hpp"
#include "geometry.hpp"
#include "material.hpp"
#include "tally.hpp"
#include <map>
#include <string>

class Simulation {
private:
    int N_particles;
    MaterialInfo mat_box;
    MaterialInfo mat_outside;
    BoxBoundary box;
    WorldBoundary world;
    int max_history_save;

    std::map<int, std::vector<double>> E_data_total;
    std::map<int, std::vector<double>> Sig_data_total;
    std::map<int, std::vector<double>> E_data_scatter;
    std::map<int, std::vector<double>> Sig_data_scatter;

    Tally results;

public:
    Simulation(int N, const std::string& mat_name, BoxBoundary b, WorldBoundary w, int max_save = 50);

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