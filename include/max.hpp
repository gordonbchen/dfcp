#pragma once

#include <vector>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"


void max_step(Clusters& clusters, const std::vector<char>& x, const HyperParams& HP, const Params& params);

void add_seqs(Clusters& clusters, const std::vector<char>& x, const Params& params, HyperParams& HP);

