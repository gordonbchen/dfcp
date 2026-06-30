#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "hyperparams.hpp"


struct Mode {
    int idx;
    size_t count;
};

struct SoftCluster {
    size_t n;
    const bool is_r;
    const int l;
    std::vector<size_t> nk;

    std::vector<SoftCluster*> parents;
    std::vector<SoftCluster*> children;

    SoftCluster(size_t n_, bool is_r_, int l_, std::vector<size_t> nk_);

    void add_child(SoftCluster *child);

    Mode mode();
};


struct SoftClusters {
    std::unordered_map<SoftCluster*, std::unique_ptr<SoftCluster>> all_clusters;
    const HyperParams& HP;
    std::vector<SoftCluster*> r_assign;
    std::vector<SoftCluster*> q_assign;
    std::vector<std::unordered_set<SoftCluster*>> rs;
    std::vector<std::unordered_set<SoftCluster*>> qs;
    int nR;

    SoftClusters(const HyperParams& HP_, const std::vector<char>& x);

    SoftCluster* create_cluster(const std::vector<int>& seqs, const std::vector<char>& x, bool is_r, int l);

    void cluster_add(SoftCluster* cluster, int idx, int emission);

    void cluster_remove(SoftCluster* cluster, int idx, int emission);
};

