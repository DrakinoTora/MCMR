#include "exporter.hpp"
#include <fstream>
#include <iostream>

void export_to_xml(const Tally& tally, const std::string& filename) {
    std::ofstream f(filename);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<mcmr_results>\n";
    f << "  <summary>\n";
    f << "    <transmission>" << tally.transmission << "</transmission>\n";
    f << "    <time_taken_seconds>" << tally.time_taken << "</time_taken_seconds>\n";
    f << "  </summary>\n";

    // Region-based tally grid: independent resolution grid over the whole world.
    // Every fission/absorption event added 1 point to the cell it occurred in.
    // One <row> per iy_t (index="0" = bottom row, y=0..dy), each a comma-separated
    // list of tally_nx ints (leftmost value = ix_t=0, x=0..dx).
    f << "  <fission_tally nx=\"" << tally.tally_nx << "\" ny=\"" << tally.tally_ny
      << "\" dx=\"" << tally.tally_dx << "\" dy=\"" << tally.tally_dy << "\">\n";
    for (int iy_t = 0; iy_t < tally.tally_ny; ++iy_t) {
        f << "    <row index=\"" << iy_t << "\">";
        for (int ix_t = 0; ix_t < tally.tally_nx; ++ix_t) {
            f << tally.fission_tally[static_cast<size_t>(iy_t) * tally.tally_nx + ix_t];
            if (ix_t + 1 != tally.tally_nx) f << ",";
        }
        f << "</row>\n";
    }
    f << "  </fission_tally>\n";

    f << "  <absorp_tally nx=\"" << tally.tally_nx << "\" ny=\"" << tally.tally_ny
      << "\" dx=\"" << tally.tally_dx << "\" dy=\"" << tally.tally_dy << "\">\n";
    for (int iy_t = 0; iy_t < tally.tally_ny; ++iy_t) {
        f << "    <row index=\"" << iy_t << "\">";
        for (int ix_t = 0; ix_t < tally.tally_nx; ++ix_t) {
            f << tally.absorp_tally[static_cast<size_t>(iy_t) * tally.tally_nx + ix_t];
            if (ix_t + 1 != tally.tally_nx) f << ",";
        }
        f << "</row>\n";
    }
    f << "  </absorp_tally>\n";
    
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

    f << "  <group_born>\n    ";
    for (size_t i = 0; i < tally.G_born.size(); ++i) {
        f << tally.G_born[i] << (i + 1 == tally.G_born.size() ? "" : ",");
    }
    f << "\n  </group_born>\n";

    f << "  <group_leak>\n    ";
    for (size_t i = 0; i < tally.G_leak.size(); ++i) {
        f << tally.G_leak[i] << (i + 1 == tally.G_leak.size() ? "" : ",");
    }
    f << "\n  </group_leak>\n";

    f << "  <trajectories>\n";
    for (size_t i = 0; i < tally.x_history.size(); ++i) {
        int gen_tag = (i < tally.history_generation.size()) ? tally.history_generation[i] : 1;
        f << "    <particle_history id=\"" << i << "\" generation=\"" << gen_tag << "\">\n";
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