#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <vector>
#include "util.hpp"


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

void parse_double(char *s, double& x) {
    char* end_ptr = nullptr;
    x = std::strtod(s, &end_ptr);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse double arg value."); };
}

void parse_int(char *s, int& x) {
    char* end_ptr = nullptr;
    x = std::strtol(s, &end_ptr, 10);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse int arg value."); };
}

