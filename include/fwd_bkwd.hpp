#pragma once
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"
#include "seq_array.hpp"


void normalize_ll(std::vector<double>& ll, int L, int K);

struct FwdBkwdMsgs {
    std::vector<double> a;
    std::vector<double> new_a;

    FwdBkwdMsgs(std::uint32_t n_cluster_ids, int L);
};

struct FwdBkwdBuffers {
    FwdBkwdMsgs bkwd;
    FwdBkwdMsgs fwd;

    FwdBkwdBuffers(std::uint32_t n_cluster_ids, int L);
};

void fwd_bkwd(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    FwdBkwdBuffers& fwd_bkwd_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);
