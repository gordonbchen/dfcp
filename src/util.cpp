#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <utility>
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

std::pair<std::vector<int8_t>, std::vector<size_t>> count_minor_alleles(
    const std::vector<size_t>& emission_counts, const std::vector<int>& masked_ls, const int n_masked_ls, const int K
) {
    std::vector<int8_t> minor_alleles(n_masked_ls);
    std::vector<size_t> minor_allele_counts(n_masked_ls);
    for (int il = 0; il < n_masked_ls; ++il) {
        auto l_counts = emission_counts.begin() + masked_ls[il]*K;
        auto min_it = std::min_element(l_counts, l_counts + K);
        minor_alleles[il] = std::distance(l_counts, min_it);
        minor_allele_counts[il] = l_counts[minor_alleles[il]];
    }
    return {std::move(minor_alleles), std::move(minor_allele_counts)};
}


double parse_double(char *s) {
    char* end_ptr = nullptr;
    double x = std::strtod(s, &end_ptr);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse double arg value."); };
    return x;
}

int parse_int(char *s) {
    char* end_ptr = nullptr;
    int x = std::strtol(s, &end_ptr, 10);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse int arg value."); };
    return x;
}

