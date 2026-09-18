#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include "clusters.hpp"
#include "model_array.hpp"
#include "seq_array.hpp"
#include "hyperparams.hpp"
#include "pbwt.hpp"
#include "util.hpp"


Cluster::Cluster(uint32_t id_, bool is_r_, int l_, int emission_)
    : id(id_), is_r(is_r_), l(l_), emission(emission_),
      n(0), q_parent(nullptr), q_child(nullptr) {}

void Cluster::add_child(Cluster *child) {
    if (is_r) {
        if (child->q_parent != nullptr) {
            throw std::runtime_error("Child q cluster has a parent already.");
        }
        children.push_back(child);
        child->q_parent = this;
        return;
    }
    if (q_child != nullptr) {
        throw std::runtime_error("Parent q cluster has a child already.");
    }
    q_child = child;
    child->parents.push_back(this);
}

Clusters::Clusters(const HyperParams& HP_) :
    next_cluster_id(0),
    HP(HP_),
    nR(0),
    rs(HP.L), qs(HP.L-1),
    rs_by_emit(HP_.total_emissions())
{}

struct PairPointerHash {
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T*, U*>& p) const {
        auto hash1 = std::hash<T*>{}(p.first);
        auto hash2 = std::hash<U*>{}(p.second);
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
};

void init_q_clusters(Clusters& clusters) {
    const HyperParams& HP = clusters.HP;
    for (int l = 0; l < HP.L-1; ++l) {
        std::unordered_map<std::pair<Cluster*, Cluster*>, Cluster*, PairPointerHash> q_map;
        for (int i = 0; i < HP.N; ++i) {
            Cluster* c = clusters.r_assign[idx2d(i, l, HP.L)];
            Cluster* c_next = clusters.r_assign[idx2d(i, l + 1, HP.L)];
            std::pair<Cluster*, Cluster*> cs{c, c_next};

            Cluster* q;
            auto existing = q_map.find(cs);
            if (existing != q_map.end()) {
                q = existing->second;
            }
            else {
                q = clusters.create_empty_cluster(false, l, -1);
                q_map.emplace(cs, q);
                c->add_child(q);
                q->add_child(c_next);
            }
            clusters.cluster_add(q, i);
        }
    }
}

void Clusters::emission_init(const ModelArray& x) {
    r_assign.resize(static_cast<std::size_t>(HP.N) * HP.L, nullptr);
    q_assign.resize(static_cast<std::size_t>(HP.N) * (HP.L - 1), nullptr);

    for (int l = 0; l < HP.L; ++l) {
        std::vector<Cluster*> r_by_emission(HP.n_emissions(l), nullptr);
        for (int i = 0; i < HP.N; ++i) {
            int emission = x(i, l);
            Cluster*& c = r_by_emission[emission];
            if (c == nullptr) {
                c = create_empty_cluster(true, l, emission);
            }
            cluster_add(c, i);
        }
    }
    init_q_clusters(*this);
}

void Clusters::pbwt_init(const SeqArray& x, int match_len) {
    r_assign.resize(HP.N * HP.L, nullptr);
    q_assign.resize(HP.N * (HP.L-1), nullptr);

    auto [a, d] = pbwt(x);
    std::vector<int> group(HP.N, 0);
    std::vector<Cluster*> r_by_group(HP.N, nullptr);

    for (int l = 0; l < HP.L; ++l) {
        int start = std::max(0, l-match_len+1);
        int end = std::min(HP.L, l+match_len);
        int group_idx = 0;

        for (int i = 0; i < HP.N; ++i) {
            if (i != 0 && d[idx2d(i,end,HP.L+1)] > start) {
                ++group_idx;
            }
            group[a[idx2d(i,end,HP.L+1)]] = group_idx;
        }

        std::fill(r_by_group.begin(), r_by_group.end(), nullptr);
        for (int i = 0; i < HP.N; ++i) {
            Cluster*& c = r_by_group[group[i]];
            if (c == nullptr) {
                c = create_empty_cluster(true, l, x(i, l));
            }
            cluster_add(c, i);
        }
    }

    init_q_clusters(*this);
}

Cluster* Clusters::create_empty_cluster(bool is_r, int l, int emission) {
    if (is_r != (emission != -1)) {
        throw std::invalid_argument("only r cluster can have emissions.");
    }

    uint32_t id;
    if (free_cluster_ids.empty()) {
        id = next_cluster_id++;
    }
    else {
        id = free_cluster_ids.back();
        free_cluster_ids.pop_back();
    }

    std::unique_ptr<Cluster> u_ptr = std::make_unique<Cluster>(id, is_r, l, emission);
    Cluster* ptr = u_ptr.get();
    if (id == all_clusters.size()) { all_clusters.push_back(std::move(u_ptr)); }
    else { all_clusters[id] = std::move(u_ptr); }

    if (!is_r) {
        qs[l].push_back(ptr);
        return ptr;
    }

    rs[l].push_back(ptr);
    ++nR;
    rs_by_emit[HP.emission_idx(l, emission)].push_back(ptr);
    return ptr;
}

void Clusters::cluster_add(Cluster* cluster, int idx) {
    ++cluster->n;

    if (cluster->is_r) {
        if (r_assign[idx2d(idx, cluster->l, HP.L)] != nullptr) {
            throw std::runtime_error("seq already assigned to r cluster");
        };
        r_assign[idx2d(idx, cluster->l, HP.L)] = cluster;
        return;
    }
    if (q_assign[idx2d(idx, cluster->l, HP.L-1)] != nullptr) {
        throw std::runtime_error("seq already assigned to q cluster");
    };
    q_assign[idx2d(idx, cluster->l, HP.L-1)] = cluster;
}

void erase_cluster(std::vector<Cluster*>& clusters, Cluster* cluster) {
    auto it = std::find(clusters.begin(), clusters.end(), cluster);
    *it = clusters.back();
    clusters.pop_back();
}

void Clusters::cluster_remove(Cluster* cluster, int idx) {
    --cluster->n;

    if (cluster->is_r) {
        r_assign[idx2d(idx, cluster->l, HP.L)] = nullptr;
    }
    else {
        q_assign[idx2d(idx, cluster->l, HP.L-1)] = nullptr;
    }

    if (cluster->n > 0) {
        return;
    }

    if (cluster->is_r) {
        for (Cluster* b : cluster->parents) {
            b->q_child = nullptr;
        }
        for (Cluster* b : cluster->children) {
            b->q_parent = nullptr;
        }
        erase_cluster(rs[cluster->l], cluster);
        --nR;
        erase_cluster(rs_by_emit[HP.emission_idx(cluster->l, cluster->emission)], cluster);
    }
    else {
        if (cluster->q_parent != nullptr) {
            erase_cluster(cluster->q_parent->children, cluster);
        }
        if (cluster->q_child != nullptr) {
            erase_cluster(cluster->q_child->parents, cluster);
        }
        erase_cluster(qs[cluster->l], cluster);
    }

    free_cluster_ids.push_back(cluster->id);
    all_clusters[cluster->id].reset();
}
