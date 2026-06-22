#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include "clusters.hpp"
#include "hyperparams.hpp"
#include "util.hpp"


Cluster::Cluster(size_t n_, bool is_r_, int l_, int emission_) : n(n_), is_r(is_r_), l(l_), emission(emission_) {
    if (is_r != (emission != -1)) { throw std::invalid_argument("only r cluster can have emissions."); };
}

void Cluster::add_child(Cluster *child) {
    children.push_back(child);
    child->parents.push_back(this);
}


Clusters::Clusters(const HyperParams& HP_, const std::vector<char>& x) :
    HP(HP_),
    r_assign(HP.N * HP.L, nullptr),
    q_assign(HP.N * (HP.L-1), nullptr),
    rs(HP.L), qs(HP.L-1),
    rs_by_emit(HP.L * HP.K),
    nR(0)
{
    // Block init.
    std::vector<int> seqs;
    seqs.reserve(HP.N);
    for (int i = 0; i < HP.N; ++i) {
        seqs.push_back(i);
    }

    std::vector<char> modes(HP.L);
    count_modes(modes, x, HP.N, HP.L, HP.K);

    Cluster* r = create_cluster(seqs, true, 0, modes[0]);
    Cluster* q = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        q = create_cluster(seqs, false, l, -1);
        r->add_child(q);

        r = create_cluster(seqs, true, l+1, modes[l+1]);
        q->add_child(r);
    }
}

Cluster* Clusters::create_cluster(const std::vector<int>& seqs, bool is_r, int l, int emission) {
    std::unique_ptr<Cluster> u_ptr = std::make_unique<Cluster>(seqs.size(), is_r, l, emission);
    Cluster* ptr = u_ptr.get();
    all_clusters[ptr] = std::move(u_ptr);

    if (!is_r) {
        for (const int& i : seqs) {
            q_assign[idx2d(i, l, HP.L-1)] = ptr;
        }
        qs[l].insert(ptr);
        return ptr;
    }

    for (const int& i : seqs) {
        r_assign[idx2d(i, l, HP.L)] = ptr;
    }
    rs[l].insert(ptr);
    rs_by_emit[idx2d(l, emission, HP.K)].insert(ptr);
    ++nR;
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

    // Delete cluster.
    for (Cluster* parent : cluster->parents) {
        std::swap(*std::find(parent->children.begin(), parent->children.end(), cluster), parent->children.back());
        parent->children.pop_back();
    }
    for (Cluster* child: cluster->children) {
        std::swap(*std::find(child->parents.begin(), child->parents.end(), cluster), child->parents.back());
        child->parents.pop_back();
    }

    if (cluster->is_r) {
        rs[cluster->l].erase(cluster);
        rs_by_emit[idx2d(cluster->l, cluster->emission, HP.K)].erase(cluster);
        --nR;
    }
    else {
        qs[cluster->l].erase(cluster);
    }

    all_clusters.erase(cluster);
}

int Clusters::cluster_mode(int l) {
    int max_k = 0;
    size_t max_nk = rs_by_emit[idx2d(l, 0, HP.K)].size();
    for (int k = 1; k < HP.K; ++k) {
        if (rs_by_emit[idx2d(l, k, HP.K)].size() > max_nk) {
            max_k = k;
            max_nk = rs_by_emit[idx2d(l, k, HP.K)].size();
        }
    }
    return max_k;
}

