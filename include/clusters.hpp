#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "hyperparams.hpp"


struct Mode {
    int idx;
    size_t count;
};


struct Cluster {
    const bool is_r;
    const int l;
    const bool soft;
    const int emission;

    size_t n;
    std::vector<size_t> nk;
    size_t n_obs;

    std::vector<Cluster*> parents;
    std::vector<Cluster*> children;


    Cluster(bool is_r_, int l_, bool soft_, int emission_, int K);

    void add_child(Cluster *child);

    Mode mode();

    int get_imputed_emission();
};


struct Clusters {
    std::unordered_map<Cluster*, std::unique_ptr<Cluster>> all_clusters;

    const HyperParams& HP;
    const bool soft;
    int nR;

    std::vector<Cluster*> r_assign;
    std::vector<Cluster*> q_assign;
    std::vector<std::unordered_set<Cluster*>> rs;
    std::vector<std::unordered_set<Cluster*>> qs;
    std::vector<std::unordered_set<Cluster*>> rs_by_emit;


    Clusters(const HyperParams& HP_, bool soft_, const std::vector<int8_t>& x);

    Cluster* create_cluster(
        const std::vector<int>& seqs, const std::vector<int8_t>& x, bool is_r, int l, int emission
    );

    Cluster* create_empty_cluster(bool is_r, int l, int emission);

    void cluster_add(Cluster* cluster, int idx, int emission);

    void cluster_remove(Cluster* cluster, int idx, int emission);

    int cluster_mode(int l);
};

