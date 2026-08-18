#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <vector>
#include "util.hpp"
#include "seq_array.hpp"


std::vector<size_t> count_emissions(const SeqArray& x, int K) {
    std::vector<size_t> counts(x.L * K, 0);
    for (int i = 0; i < x.N; ++i) {
        for (int l = 0; l < x.L; ++l) {
            ++counts[idx2d(l, x(i, l), K)];
        }
    }
    return counts;
}

std::vector<int8_t> get_emission_modes(const std::vector<size_t>& emission_counts, const int L, const int K) {
    std::vector<int8_t> modes(L);
    for (int l = 0; l < L; ++l) {
        auto l_counts = emission_counts.begin() + l*K;
        auto max_it = std::max_element(l_counts, l_counts + K);
        if (*max_it == 0) { throw std::runtime_error("No valid alleles at loc."); };
        modes[l] = std::distance(l_counts, max_it);
    }
    return modes;
}
