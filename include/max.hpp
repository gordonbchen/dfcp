#pragma once

#include <cstdint>
#include <vector>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"


void max_step(Clusters& clusters, const std::vector<int8_t>& x, const HyperParams& HP, const Params& params);

void add_seqs(
    Clusters& clusters, std::vector<int8_t>::const_iterator new_x, int n_new,
    HyperParams& HP, const Params& params
);

void max_cluster_emissions(Clusters& clusters, const HyperParams& HP, const Params& params);

