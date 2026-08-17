#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>


inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

inline int8_t get_xil(
    std::vector<int8_t>::const_iterator xi, int l, const std::unordered_map<int, int> *unmasked_ls
) {
    if (unmasked_ls == nullptr) {
        return xi[l];
    }
    return unmasked_ls->contains(l) ? xi[unmasked_ls->at(l)] : -1;
}

std::vector<size_t> count_emissions(const std::vector<int8_t>& x, const int N, const int L, const int K);
std::vector<int8_t> get_emission_modes(const std::vector<size_t>& emission_counts, const int L, const int K);

