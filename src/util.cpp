#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>
#include "util.hpp"


size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

std::vector<char> count_modes(const std::vector<char>& x, const int N, const int L, const int K) {
    std::vector<char> modes(L);
    std::vector<int> counts(K);
    for (int l = 0; l < L; ++l) {
        std::fill(counts.begin(), counts.end(), 0);
        for (int i = 0; i < N; ++i) {
            if (x[idx2d(i, l, L)] != -1) {
                ++counts[x[idx2d(i, l, L)]];
            }
        }
        auto max_it = std::max_element(counts.begin(), counts.end());
        if (*max_it == 0) { throw std::runtime_error("No valid alleles at loc."); };
        modes[l] = std::distance(counts.begin(), max_it);
    }
    return modes;
}

