#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include "soft_clusters.hpp"
#include "hyperparams.hpp"
#include "util.hpp"


SoftCluster::SoftCluster(size_t n_, bool is_r_, int l_, std::vector<size_t> nk_) : n(n_), is_r(is_r_), l(l_), nk(nk_) {
    if (is_r != (nk.size() > 0 )) { throw std::invalid_argument("only r cluster can have nk."); };
}

void SoftCluster::add_child(SoftCluster *child) {
    children.push_back(child);
    child->parents.push_back(this);
}

Mode SoftCluster::mode() {
    Mode mode{0, nk[0]};
    for (size_t k = 1; k < nk.size(); ++k) {
        if (nk[k] > mode.count) {
            mode.idx = k;
            mode.count = nk[k];
        }
    }
    return mode;
}


SoftClusters::SoftClusters(const HyperParams& HP_, const std::vector<char>& x) :
    HP(HP_),
    r_assign(HP.N * HP.L, nullptr),
    q_assign(HP.N * (HP.L-1), nullptr),
    rs(HP.L), qs(HP.L-1),
    nR(0)
{
    // Block init.
    std::vector<int> seqs;
    seqs.reserve(HP.N);
    for (int i = 0; i < HP.N; ++i) {
        seqs.push_back(i);
    }

    SoftCluster* r = create_cluster(seqs, x, true, 0);
    SoftCluster* q = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        q = create_cluster(seqs, x, false, l);
        r->add_child(q);

        r = create_cluster(seqs, x, true, l+1);
        q->add_child(r);
    }
}

std::vector<size_t> count_emissions(const std::vector<int>& seqs, const std::vector<char>& x, int l, int L, int K) {
    std::vector<size_t> nk(K, 0);
    for (int i : seqs) {
        if (x[idx2d(i, l, L)] != -1) {
            ++nk[x[idx2d(i, l, L)]];
        }
    }
    return nk;
}

SoftCluster* SoftClusters::create_cluster(const std::vector<int>& seqs, const std::vector<char>& x, bool is_r, int l) {
    std::vector<size_t> nk = is_r ? count_emissions(seqs, x, l, HP.L, HP.K) : std::vector<size_t>(0);
    std::unique_ptr<SoftCluster> u_ptr = std::make_unique<SoftCluster>(seqs.size(), is_r, l, std::move(nk));
    SoftCluster* ptr = u_ptr.get();
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
    return ptr;
}

void SoftClusters::cluster_add(SoftCluster* cluster, int idx, int emission) {
    ++cluster->n;

    if (cluster->is_r) {
        if (r_assign[idx2d(idx, cluster->l, HP.L)] != nullptr) {
            throw std::runtime_error("seq already assigned to r cluster");
        };
        r_assign[idx2d(idx, cluster->l, HP.L)] = cluster;
        if (emission != -1) {
            ++cluster->nk[emission];
        }
        return;
    }
    if (q_assign[idx2d(idx, cluster->l, HP.L-1)] != nullptr) {
        throw std::runtime_error("seq already assigned to q cluster");
    };
    q_assign[idx2d(idx, cluster->l, HP.L-1)] = cluster;
}

void SoftClusters::cluster_remove(SoftCluster* cluster, int idx, int emission) {
    --cluster->n;

    if (cluster->is_r) {
        r_assign[idx2d(idx, cluster->l, HP.L)] = nullptr;
        if (emission != -1) {
            --cluster->nk[emission];
        }
    }
    else {
        q_assign[idx2d(idx, cluster->l, HP.L-1)] = nullptr;
    }

    if (cluster->n > 0) {
        return;
    }

    // Delete cluster.
    for (SoftCluster* parent : cluster->parents) {
        std::swap(*std::find(parent->children.begin(), parent->children.end(), cluster), parent->children.back());
        parent->children.pop_back();
    }
    for (SoftCluster* child: cluster->children) {
        std::swap(*std::find(child->parents.begin(), child->parents.end(), cluster), child->parents.back());
        child->parents.pop_back();
    }

    if (cluster->is_r) {
        rs[cluster->l].erase(cluster);
        --nR;
    }
    else {
        qs[cluster->l].erase(cluster);
    }

    all_clusters.erase(cluster);
}

