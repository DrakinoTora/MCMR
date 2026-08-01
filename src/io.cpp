#include "io.hpp"
#include <fstream>
#include <iostream>

void simpan_hasil_xml(const SimulationResult& res, const std::string& filename) {
    std::ofstream file(filename);
    file << "<simulasi_mc>\n";
    file << "  <absorpsi_A>" << res.abs_A << "</absorpsi_A>\n";
    file << "  <absorpsi_B>" << res.abs_B << "</absorpsi_B>\n";
    file << "  <transmisi>" << res.transmisi << "</transmisi>\n";
    file << "  <waktu_detik>" << res.time_taken << "</waktu_detik>\n";
    file << "</simulasi_mc>\n";
    file.close();
    std::cout << "Hasil berhasil disimpan ke " << filename << std::endl;
}