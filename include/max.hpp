#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"


double get_cluster_emission_ll(
    Cluster* a, int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

std::vector<Cluster*> get_viterbi_clusters(
    std::vector<int8_t>::const_iterator xi, const std::unordered_map<int, int> *unmasked_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);


int get_new_cluster_emission(
    int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

void max_step(
    const std::vector<int8_t>& x,
    Clusters& clusters, const Params& params, const HyperParams& HP
);

void add_seqs(
    std::vector<int8_t>::const_iterator new_x, int n_new,
    Clusters& clusters, const Params& params, HyperParams& HP
);


void max_cluster_emissions(Clusters& clusters, const Params& params, const HyperParams& HP);

