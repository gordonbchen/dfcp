#pragma once

#include <unordered_map>
#include <vector>

#include "seq_array.hpp"


struct RefFilter {
    const SeqArray& ref;
    std::vector<int> observed_loci;
    std::vector<double> weight_zero;
    std::vector<double> weight_one;

    RefFilter(const SeqArray& ref, const std::unordered_map<int, int>& obs_ls);

    std::vector<int> top_references(const SeqArray& target, int target_idx, int top_k) const;
};

SeqArray select_references(const SeqArray& ref, const std::vector<int>& indexes);
