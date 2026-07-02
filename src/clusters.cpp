#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include "clusters.hpp"
#include "hyperparams.hpp"
#include "util.hpp"


Cluster::Cluster(
    bool is_r_, int l_, bool soft_, size_t n_,
    int emission_, std::vector<size_t> nk_, size_t n_obs_
)
    : is_r(is_r_), l(l_), soft(soft_), n(n_), emission(emission_), nk(nk_), n_obs(n_obs_) {}

void Cluster::add_child(Cluster *child) {
    children.push_back(child);
    child->parents.push_back(this);
}

Mode Cluster::mode() {
    Mode mode{0, nk[0]};
    for (size_t k = 1; k < nk.size(); ++k) {
        if (nk[k] > mode.count) {
            mode.idx = k;
            mode.count = nk[k];
        }
    }
    return mode;
}

int Cluster::get_imputed_emission() {
    return soft ? mode().idx : emission;
}


Clusters::Clusters(const HyperParams& HP_, bool soft_, const std::vector<int8_t>& x) :
    HP(HP_),
    soft(soft_),
    nR(0),
    r_assign(HP.N * HP.L, nullptr), q_assign(HP.N * (HP.L-1), nullptr),
    rs(HP.L), qs(HP.L-1),
    rs_by_emit(soft_ ? 0 : HP_.L * HP_.K)
{
    // Block init.
    std::vector<int> seqs;
    seqs.reserve(HP.N);
    for (int i = 0; i < HP.N; ++i) {
        seqs.push_back(i);
    }

    std::vector<int8_t> modes{count_modes(x, HP.N, HP.L, HP.K)};

    Cluster* r = create_cluster(seqs, x, true, 0, modes[0]);
    Cluster* q = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        q = create_cluster(seqs, x, false, l, -1);
        r->add_child(q);

        r = create_cluster(seqs, x, true, l+1, modes[l+1]);
        q->add_child(r);
    }
}

std::tuple<std::vector<size_t>, size_t> count_emissions(
    const std::vector<int>& seqs, const std::vector<int8_t>& x, const HyperParams& HP, int l
) {
    size_t n_obs = 0;
    std::vector<size_t> nk(HP.K, 0);
    for (int i : seqs) {
        if (x[idx2d(i, l, HP.L)] != -1) {
            ++nk[x[idx2d(i, l, HP.L)]];
            ++n_obs;
        }
    }
    return std::tuple<std::vector<size_t>, size_t>{nk, n_obs};
}

Cluster* Clusters::create_cluster(
    const std::vector<int>& seqs, const std::vector<int8_t>& x, bool is_r, int l, int emission
) {
    if (!soft && (is_r != (emission != -1))) {
        throw std::invalid_argument("only r cluster can have emissions.");
    }

    auto [nk, n_obs] = (soft && is_r) ? count_emissions(seqs, x, HP, l)
        : std::tuple<std::vector<size_t>, size_t>{std::vector<size_t>(0), 0};
    std::unique_ptr<Cluster> u_ptr = std::make_unique<Cluster>(
        is_r, l, soft, seqs.size(), emission, std::move(nk), n_obs
    );

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
    ++nR;
    if (!soft) {
        rs_by_emit[idx2d(l, emission, HP.K)].insert(ptr);
    }
    return ptr;
}

void Clusters::cluster_add(Cluster* cluster, int idx, int emission) {
    ++cluster->n;

    if (cluster->is_r) {
        if (r_assign[idx2d(idx, cluster->l, HP.L)] != nullptr) {
            throw std::runtime_error("seq already assigned to r cluster");
        };
        r_assign[idx2d(idx, cluster->l, HP.L)] = cluster;
        if (soft && (emission != -1)) {
            ++cluster->nk[emission];
            ++cluster->n_obs;
        }
        return;
    }
    if (q_assign[idx2d(idx, cluster->l, HP.L-1)] != nullptr) {
        throw std::runtime_error("seq already assigned to q cluster");
    };
    q_assign[idx2d(idx, cluster->l, HP.L-1)] = cluster;
}

void Clusters::cluster_remove(Cluster* cluster, int idx, int emission) {
    --cluster->n;

    if (cluster->is_r) {
        r_assign[idx2d(idx, cluster->l, HP.L)] = nullptr;
        if (soft && (emission != -1)) {
            --cluster->nk[emission];
            --cluster->n_obs;
        }
    }
    else {
        q_assign[idx2d(idx, cluster->l, HP.L-1)] = nullptr;
    }

    if (cluster->n > 0) {
        return;
    }

    // Delete cluster.
    for (Cluster* parent : cluster->parents) {
        std::swap(
            *std::find(parent->children.begin(), parent->children.end(), cluster), parent->children.back()
        );
        parent->children.pop_back();
    }
    for (Cluster* child: cluster->children) {
        std::swap(*std::find(child->parents.begin(), child->parents.end(), cluster), child->parents.back());
        child->parents.pop_back();
    }

    if (cluster->is_r) {
        rs[cluster->l].erase(cluster);
        --nR;
        if (!soft) {
            rs_by_emit[idx2d(cluster->l, cluster->emission, HP.K)].erase(cluster);
        }
    }
    else {
        qs[cluster->l].erase(cluster);
    }

    all_clusters.erase(cluster);
}

int Clusters::cluster_mode(int l) {
    if (soft) { throw std::runtime_error("Soft clusters don't have emissions, no cluster emission mode"); }
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

