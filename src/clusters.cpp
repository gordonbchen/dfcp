#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include "clusters.hpp"
#include "hyperparams.hpp"
#include "util.hpp"


Cluster::Cluster(bool is_r_, int l_, int emission_, int K)
    : is_r(is_r_), l(l_), emission(emission_),
      n(0), nk(K, 0), n_obs(0) {}

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

int Cluster::get_imputed_emission(bool soft) {
    return soft ? mode().idx : emission;
}


Clusters::Clusters(const HyperParams& HP_, bool soft_, bool noisy_) :
    HP(HP_),
    soft(soft_), noisy(noisy_),
    nR(0),
    rs(HP.L), qs(HP.L-1),
    rs_by_emit(soft_ ? 0 : HP_.L * HP_.K),
    n_matches(0), n_obs(0)
{}

void Clusters::block_init(const std::vector<int8_t>& x) {
    r_assign.resize(HP.N * HP.L, nullptr);
    q_assign.resize(HP.N * (HP.L-1), nullptr);

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

Cluster* Clusters::create_cluster(
    const std::vector<int>& seqs, const std::vector<int8_t>& x, bool is_r, int l, int emission
) {
    Cluster* c = create_empty_cluster(is_r, l, emission);
    for (const int& i : seqs) {
        cluster_add(c, i, x[idx2d(i, l, HP.L)]);
    }
    return c;
}

Cluster* Clusters::create_empty_cluster(bool is_r, int l, int emission) {
    if (!soft && (is_r != (emission != -1))) {
        throw std::invalid_argument("only r cluster can have emissions.");
    }

    std::unique_ptr<Cluster> u_ptr = std::make_unique<Cluster>(is_r, l, emission, HP.K);
    Cluster* ptr = u_ptr.get();
    all_clusters[ptr] = std::move(u_ptr);

    if (!is_r) {
        qs[l].insert(ptr);
        return ptr;
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
        if ((soft || noisy) && (emission != -1)) {
            ++cluster->nk[emission];
            ++cluster->n_obs;
        }
        if (noisy && (emission != -1)) {
            ++n_obs;
            if (emission == cluster->emission) {
                ++n_matches;
            }
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
        if ((soft || noisy) && (emission != -1)) {
            --cluster->nk[emission];
            --cluster->n_obs;
        }
        if (noisy && (emission != -1)) {
            --n_obs;
            if (emission == cluster->emission) {
                --n_matches;
            }
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

void Clusters::set_emission(Cluster* c, int new_emission) {
    if (!noisy || !c->is_r) {
        throw std::runtime_error("set_emission can only be called for noisy R clusters.");
    }
    if (c->emission == new_emission) {
        return;
    }
    rs_by_emit[idx2d(c->l,c->emission,HP.K)].erase(c);
    rs_by_emit[idx2d(c->l,new_emission,HP.K)].insert(c);
    n_matches += c->nk[new_emission];
    n_matches -= c->nk[c->emission];
    c->emission = new_emission;
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

