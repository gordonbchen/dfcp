#pragma once
#include <vector>
#include <cstdint>
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"


std::vector<double> forward_backward(
    std::vector<int8_t>::const_iterator xi, std::vector<int> prob_idxs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

