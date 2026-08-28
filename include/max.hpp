#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "seq_array.hpp"


double get_cluster_emission_ll(
    Cluster* a, int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);


struct ViterbiMsg {
    double ll;
    Cluster* next;
};

std::vector<Cluster*> get_viterbi_clusters(
    const SeqArray& x, int i, const std::unordered_map<int, int> *obs_ls,
    std::vector<std::unordered_map<Cluster*, ViterbiMsg>>& a_msgs,
    std::vector<std::unordered_map<Cluster*, ViterbiMsg>>& b_msgs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);


int get_new_cluster_emission(
    int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

void max_step(const SeqArray& x, Clusters& clusters, const Params& params, const HyperParams& HP);

void add_seqs(const SeqArray& x, Clusters& clusters, const Params& params, HyperParams& HP);


void max_cluster_emissions(Clusters& clusters, const Params& params, const HyperParams& HP);
