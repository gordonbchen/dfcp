#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>
#include "hyperparams.hpp"
#include "model_array.hpp"
#include "seq_array.hpp"


struct Cluster {
    const uint32_t id;
    const bool is_r;
    const int l;
    const int emission;

    size_t n;

    std::vector<Cluster*> parents;
    std::vector<Cluster*> children;
    Cluster* q_parent;
    Cluster* q_child;


    Cluster(uint32_t id_, bool is_r_, int l_, int emission_);

    void add_child(Cluster *child);
};


struct Clusters {
    std::vector<std::unique_ptr<Cluster>> all_clusters;
    uint32_t next_cluster_id;
    std::vector<uint32_t> free_cluster_ids;

    const HyperParams& HP;
    int nR;

    std::vector<Cluster*> r_assign;
    std::vector<Cluster*> q_assign;
    std::vector<std::vector<Cluster*>> rs;
    std::vector<std::vector<Cluster*>> qs;

    std::vector<std::vector<Cluster*>> rs_by_emit;


    Clusters(const HyperParams& HP_);

    void emission_init(const ModelArray& x);
    void pbwt_init(const SeqArray& x, int match_len);

    Cluster* create_empty_cluster(bool is_r, int l, int emission);

    void cluster_add(Cluster* cluster, int idx);

    void cluster_remove(Cluster* cluster, int idx);

    const std::vector<Cluster*>& get_matching_as(int l, int emission) const {
        return emission == -1 ? rs[l] : rs_by_emit[HP.emission_idx(l, emission)];
    }
};
