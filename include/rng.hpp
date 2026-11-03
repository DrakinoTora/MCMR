#pragma once
#include <cstdint>
#include <optional>
#include <random>

// Seeding helpers shared by both engines.
//
// A run's randomness comes from ONE std::mt19937 built from a 64-bit seed.
// std::mt19937's own constructor only takes 32 bits (it would silently
// truncate), so the seed is spread over two 32-bit words through std::seed_seq
// -- seeds that differ only in their upper 32 bits therefore still give
// different streams.
inline std::mt19937 make_rng(std::uint64_t seed) {
    std::seed_seq ss{static_cast<std::uint32_t>(seed & 0xFFFFFFFFu),
                     static_cast<std::uint32_t>(seed >> 32)};
    return std::mt19937(ss);
}

// seed given -> use it as-is (reproducible run).
// seed empty -> draw a fresh random one. Either way the caller records the
// value actually used (Tally::seed), so ANY run can be repeated afterwards.
inline std::uint64_t resolve_seed(const std::optional<std::uint64_t>& seed) {
    if (seed) return *seed;
    std::random_device rd;
    return (static_cast<std::uint64_t>(rd()) << 32) | static_cast<std::uint64_t>(rd());
}
