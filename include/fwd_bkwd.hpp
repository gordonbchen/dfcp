#pragma once
#include <unordered_map>
#include <vector>
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"
#include "seq_array.hpp"


void normalize_ll(std::vector<double>& ll, int L, int K);

std::vector<double> fwd_bkwd(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);
