#pragma once

#include <cstdint>
#include <vector>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"


double get_cluster_emission_ll(
    Cluster* a, int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

std::vector<Cluster*> get_viterbi_clusters(
    const Clusters& clusters, std::vector<int8_t>::const_iterator xi,
    const HyperParams& HP, const Params& params
);


int get_new_cluster_emission(int8_t xil, int l, const Clusters& clusters, const HyperParams& HP, const Params& params);

void max_step(Clusters& clusters, const std::vector<int8_t>& x, const HyperParams& HP, const Params& params);

void add_seqs(
    Clusters& clusters, std::vector<int8_t>::const_iterator new_x, int n_new,
    HyperParams& HP, const Params& params
);


void max_cluster_emissions(Clusters& clusters, const HyperParams& HP, const Params& params);

