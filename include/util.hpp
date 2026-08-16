#pragma once

#include <cstddef>
#include <vector>
#include <cstdint>


inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

std::vector<size_t> count_emissions(const std::vector<int8_t>& x, const int N, const int L, const int K);
std::vector<int8_t> get_emission_modes(const std::vector<size_t>& emission_counts, const int L, const int K);

