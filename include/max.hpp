#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "clusters.hpp"
#include "hyperparams.hpp"
#include "model_array.hpp"
#include "params.hpp"
#include "seq_array.hpp"

struct ScalarObs;
struct BlockObs;


struct ViterbiMsg {
    double ll;
    Cluster* next;
};

struct ViterbiBuffers {
    std::vector<ViterbiMsg> a_msgs;
    std::vector<ViterbiMsg> new_a_msgs;
    std::vector<ViterbiMsg> new_b_msgs;
    std::vector<Cluster*> path;

    ViterbiBuffers(uint32_t n_cluster_ids, int L);
};

void get_viterbi_path(
    const ScalarObs& obs, ViterbiBuffers& viterbi_bufs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

void get_viterbi_path(
    const BlockObs& obs, ViterbiBuffers& viterbi_bufs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);


void max_step(
    const ModelArray& x, Clusters& clusters, const Params& params, const HyperParams& HP,
    int batch_size
);

void add_seqs(const ModelArray& x, Clusters& clusters, const Params& params, HyperParams& HP);


void get_viterbi_impute_probs(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    ViterbiBuffers& viterbi_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);

void get_blocked_viterbi_impute_probs(
    const ModelArray& ref, const BlockObs& obs, const std::unordered_map<int, int>& obs_ls,
    ViterbiBuffers& viterbi_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);
