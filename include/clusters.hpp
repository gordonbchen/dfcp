#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "hyperparams.hpp"


struct Cluster {
    size_t n;
    const bool is_r;
    const int l;
    const int emission;

    std::vector<Cluster*> parents;
    std::vector<Cluster*> children;

    Cluster(size_t n_, bool is_r_, int l_, int emission_);

    void add_child(Cluster *child);
};


struct Clusters {
    std::unordered_map<Cluster*, std::unique_ptr<Cluster>> all_clusters;
    const HyperParams& HP;
    std::vector<Cluster*> r_assign;
    std::vector<Cluster*> q_assign;
    std::vector<std::unordered_set<Cluster*>> rs;
    std::vector<std::unordered_set<Cluster*>> qs;
    std::vector<std::unordered_set<Cluster*>> rs_by_emit;
    int nR;

    Clusters(const HyperParams& HP_, const std::vector<char>& x);

    Cluster* create_cluster(const std::vector<int>& seqs, bool is_r, int l, int emission);

    void cluster_add(Cluster* cluster, int idx);

    void cluster_remove(Cluster* cluster, int idx);

    int cluster_mode(int l);
};

