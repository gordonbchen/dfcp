#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <vector>
#include "util.hpp"


std::vector<size_t> count_emissions(const std::vector<int8_t>& x, const int N, const int L, const int K) {
    std::vector<size_t> counts(L*K, 0);
    for (int i = 0; i < N; ++i) {
        for (int l = 0; l < L; ++l) {
            int emission = x[idx2d(i, l, L)];
            if (emission != -1) {
                ++counts[idx2d(l, emission, K)];
            }
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

