#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "hyperparams.hpp"
#include "seq_array.hpp"


struct Mode {
    int idx;
    size_t count;
};

struct Cluster {
    const bool is_r;
    const int l;
    int emission;

    size_t n;
    std::vector<size_t> nk;
    size_t n_obs;

    std::vector<Cluster*> parents;
    std::vector<Cluster*> children;


    Cluster(bool is_r_, int l_, int emission_, int K);

    void add_child(Cluster *child);

    Mode mode() const;
};


enum class EmitMode { hard, noisy, soft };

struct Clusters {
    std::unordered_map<Cluster*, std::unique_ptr<Cluster>> all_clusters;

    const HyperParams& HP;
    const EmitMode emit_mode;
    int nR;

    std::vector<Cluster*> r_assign;
    std::vector<Cluster*> q_assign;
    std::vector<std::unordered_set<Cluster*>> rs;
    std::vector<std::unordered_set<Cluster*>> qs;

    std::vector<std::unordered_set<Cluster*>> rs_by_emit;
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

    int cluster_mode(int l) const;
};
