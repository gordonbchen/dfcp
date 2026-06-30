#pragma once

#include <vector>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "soft_clusters.hpp"


template <typename TCluster>
struct Msg {
    double ll;
    TCluster* next;
};


void max_step(Clusters& clusters, const std::vector<char>& x, const HyperParams& HP, const Params& params);

void max_step(SoftClusters& clusters, const std::vector<char>& x, const HyperParams& HP, const Params& params);


void viterbi_seq(Clusters& clusters, const std::vector<char>& x, int i, const HyperParams& HP, const Params& params);

void viterbi_seq(SoftClusters& clusters, const std::vector<char>& x, int i, const HyperParams& HP, const Params& params);


template <typename TClusters>
void add_seqs(TClusters& clusters, const std::vector<char>& x, const Params& params, HyperParams& HP) {
    int n = x.size() / HP.L;
    int old_N = HP.N;
    HP.N += n;
    clusters.r_assign.resize(clusters.r_assign.size() + x.size(), nullptr);
    clusters.q_assign.resize(clusters.q_assign.size() + n * (HP.L-1), nullptr);

    for (int i = 0; i < n; ++i) {
        viterbi_seq(clusters, x, old_N + i, HP, params);
    }
}

