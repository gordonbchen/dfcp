#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <vector>
#include "util.hpp"
#include "seq_array.hpp"


double parse_double(const char* value) {
    char* end = nullptr;
    double parsed = std::strtod(value, &end);
    if (end == value) { throw std::invalid_argument("Failed to parse double arg value."); }
    return parsed;
}

int parse_int(const char* value) {
    char* end = nullptr;
    int parsed = std::strtol(value, &end, 10);
    if (end == value) { throw std::invalid_argument("Failed to parse int arg value."); }
    return parsed;
}

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
        modes[l] = std::distance(l_counts, max_it);
    }
    return modes;
}
