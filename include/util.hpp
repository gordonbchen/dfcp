#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "seq_array.hpp"


inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

inline int8_t get_xil(const SeqArray& x, int i, int l, const std::unordered_map<int, int> *unmasked_ls) {
    if (unmasked_ls == nullptr) {
        return x(i, l);
    }
    auto mapped_l = unmasked_ls->find(l);
    return mapped_l == unmasked_ls->end() ? -1 : x(i, mapped_l->second);
}

std::vector<size_t> count_emissions(const SeqArray& x, int K);
std::vector<int8_t> get_emission_modes(const std::vector<size_t>& emission_counts, const int L, const int K);
