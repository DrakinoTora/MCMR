#include "exporter.hpp"
#include <fstream>
#include <iostream>

void export_to_xml(const Tally& tally, const std::string& filename) {
    std::ofstream f(filename);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<mcmr_results>\n";
    f << "  <summary>\n";
    for (const auto& kv : tally.absorp_by_material) {
        f << "    <absorp material=\"" << kv.first << "\">" << kv.second << "</absorp>\n";
    }
    f << "    <transmisi>" << tally.transmisi << "</transmisi>\n";
    f << "    <time_taken_seconds>" << tally.time_taken << "</time_taken_seconds>\n";
    f << "  </summary>\n";
    
    f << "  <energy_born>\n    ";
    for (size_t i = 0; i < tally.E_born.size(); ++i) {
        f << tally.E_born[i] << (i + 1 == tally.E_born.size() ? "" : ",");
    }
    f << "\n  </energy_born>\n";

    f << "  <energy_leak>\n    ";
    for (size_t i = 0; i < tally.E_leak.size(); ++i) {
        f << tally.E_leak[i] << (i + 1 == tally.E_leak.size() ? "" : ",");
    }
    f << "\n  </energy_leak>\n";

    f << "  <trajectories>\n";
    for (size_t i = 0; i < tally.x_history.size(); ++i) {
        f << "    <particle_history id=\"" << i << "\">\n";
        f << "      <x>";
        for (size_t j = 0; j < tally.x_history[i].size(); ++j) {
            f << tally.x_history[i][j] << (j + 1 == tally.x_history[i].size() ? "" : ",");
        }
        f << "</x>\n      <y>";
        for (size_t j = 0; j < tally.y_history[i].size(); ++j) {
            f << tally.y_history[i][j] << (j + 1 == tally.y_history[i].size() ? "" : ",");
        }
        f << "</y>\n    </particle_history>\n";
    }
    f << "  </trajectories>\n";
    f << "</mcmr_results>\n";
    f.close();
}