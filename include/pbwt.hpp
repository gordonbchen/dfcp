#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <utility>
#include "seq_array.hpp"

struct PbwtBlocks {
    std::vector<std::uint32_t> labels;
    std::vector<std::size_t> rep_offsets;
    std::vector<int> reps;
    std::vector<std::uint32_t> emission_sizes;
    std::vector<int> snp_offsets;
    int L;
    int n_max_k_cuts;
    int n_compat_cuts;

    std::uint32_t operator()(int i, int l) const {
        return labels[static_cast<std::size_t>(i) * L + l];
    }

    int n_emissions(int l) const {
        return static_cast<int>(rep_offsets[l + 1] - rep_offsets[l]);
    }

    int rep(int l, int emission) const {
        return reps[rep_offsets[l] + emission];
    }

    std::uint32_t emission_size(int l, int emission) const {
        return emission_sizes[rep_offsets[l] + emission];
    }

    int snp_start(int l) const {
        return snp_offsets[l];
    }

    int snp_end(int l) const {
        return snp_offsets[l + 1];
    }
};


std::pair<std::vector<int>, std::vector<int>> pbwt(const SeqArray& x);
PbwtBlocks pbwt_blocks_greedy(
    const SeqArray& ref, int max_k,
    const SeqArray* target,
    const std::unordered_map<int, int>* obs_ls
);
