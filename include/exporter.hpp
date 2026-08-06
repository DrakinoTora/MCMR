#pragma once
#include "tally.hpp"
#include <string>

void export_to_xml(const Tally& tally, const std::string& filename = "mcmr_results.xml");