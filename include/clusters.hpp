#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>
#include "hyperparams.hpp"
#include "seq_array.hpp"
#include "util.hpp"


struct Mode {
    int idx;
    size_t count;
};

struct Cluster {
    const uint32_t id;
    const bool is_r;
    const int l;
    int emission;

    size_t n;
    std::vector<size_t> nk;
    size_t n_obs;

    std::vector<Cluster*> parents;
    std::vector<Cluster*> children;
    Cluster* q_parent;
    Cluster* q_child;


    Cluster(uint32_t id_, bool is_r_, int l_, int emission_, int K);

    void add_child(Cluster *child);

    Mode mode() const;
};


enum class EmitMode { hard, noisy, soft };

struct Clusters {
    std::vector<std::unique_ptr<Cluster>> all_clusters;
    uint32_t next_cluster_id;
    std::vector<uint32_t> free_cluster_ids;

    const HyperParams& HP;
    const EmitMode emit_mode;
    int nR;

    std::vector<Cluster*> r_assign;
    std::vector<Cluster*> q_assign;
    std::vector<std::vector<Cluster*>> rs;
    std::vector<std::vector<Cluster*>> qs;

    std::vector<std::vector<Cluster*>> rs_by_emit;
    size_t n_matches;
    size_t n_obs;


    Clusters(const HyperParams& HP_, EmitMode emit_mode_);

    void block_init(const SeqArray& x);
    void pbwt_init(const SeqArray& x, int match_len, bool match_curr);

    Cluster* create_cluster(const std::vector<int>& seqs, const SeqArray& x, bool is_r, int l, int emission);

    Cluster* create_empty_cluster(bool is_r, int l, int emission);

    void cluster_add(Cluster* cluster, int idx, int emission);

    void cluster_remove(Cluster* cluster, int idx, int emission);

    void set_emission(Cluster* c, int new_emission);

    const std::vector<Cluster*>& get_matching_as(int l, int emission) const {
        return emit_mode != EmitMode::hard || emission == -1 ?
            rs[l] : rs_by_emit[idx2d(l, emission, HP.K)];
    }

    int cluster_mode(int l) const;
};
