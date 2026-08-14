#pragma once
#include <vector>
#include <cstdint>
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"


void normalize_ll(std::vector<double>& ll, int L, int K);

std::vector<double> fwd_bkwd(
    std::vector<int8_t>::const_iterator xi, const std::vector<int>& prob_idxs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

