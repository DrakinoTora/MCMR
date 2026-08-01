#pragma once
#include "particle.hpp"
#include <string>

// Fungsi pendukung XML jika pengguna ingin menyimpan hasil ke disk
void simpan_hasil_xml(const SimulationResult& res, const std::string& filename = "hasil_simulasi.xml");