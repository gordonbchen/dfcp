#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include "clusters.hpp"
#include "hyperparams.hpp"
#include "params.hpp"
#include "math.hpp"
#include "util.hpp"


Cluster::Cluster(size_t n_, bool is_r_, int l_, int emission_) : n(n_), is_r(is_r_), l(l_), emission(emission_) {
    if (is_r != (emission != -1)) { throw std::invalid_argument("only r cluster can have emissions."); };
}

void Cluster::add_child(Cluster *child) {
    children.push_back(child);
    child->parents.push_back(this);
}


void count_modes(std::vector<char>& modes, const std::vector<char>& x, const int N, const int L, const int K) {
    std::vector<int> counts(K);
    for (int l = 0; l < L; ++l) {
        std::fill(counts.begin(), counts.end(), 0);
        for (int i = 0; i < N; ++i) {
            if (x[idx2d(i, l, L)] != -1) {
                ++counts[x[idx2d(i, l, L)]];
            }
        }
        auto max_it = std::max_element(counts.begin(), counts.end());
        if (*max_it == 0) { throw std::runtime_error("No valid alleles at loc."); };
        modes[l] = std::distance(counts.begin(), max_it);
    }
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

void Clusters::max_step(const std::vector<char>& x, const Params& params) {
    for (int i = 0; i < HP.N; ++i) {
        for (int l = 0; l < HP.L; ++l) {
            cluster_remove(r_assign[idx2d(i, l, HP.L)], i);
            if (l == HP.L-1) {
                break;
            }
            cluster_remove(q_assign[idx2d(i, l, HP.L-1)], i);
        }
        viterbi_seq(x.begin() + idx2d(i, 0, HP.L), i, params);
    }
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

struct Msg {
    double ll;
    Cluster* next;
};

double get_msg_ll(const std::unordered_map<Cluster*, Msg>& msgs, Cluster* c) {
    auto it = msgs.find(c);
    if (it == msgs.end()) {
        return -std::numeric_limits<double>::infinity();
    }
    return it->second.ll;
}

void Clusters::viterbi_seq(std::vector<char>::const_iterator xi, int i, const Params& params) {
    std::vector<std::unordered_map<Cluster*, Msg>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, Msg>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        // Likelihood for new cluster.
        auto& ma = a_msgs[l];
        int emission = xi[l] == -1 ? cluster_mode(l) : xi[l];
        double new_a_ll = (
            delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, rs_by_emit[idx2d(l,emission,HP.K)].size())
            - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, rs[l].size())
        );

        std::unordered_set<Cluster*>& matching_as = (
            xi[l] == -1 ? rs[l] : rs_by_emit[idx2d(l, emission, HP.K)]
        );
        if (l == HP.L-1) {
            ma[nullptr] = Msg{new_a_ll, nullptr};
            for (Cluster *a : matching_as) {
                ma[a] = Msg{0.0, nullptr};
            }
            continue;
        }

        // b messages.
        auto& next_ma = a_msgs[l+1];
        auto& mb = b_msgs[l];
        for (Cluster* b : qs[l]) {
            if (b->children.size() != 1) { throw std::runtime_error("b clusters should only have 1 child."); };
            Cluster* next_a = *b->children.begin();
            mb[b] = Msg{get_msg_ll(next_ma, next_a), next_a};
         }

        int nQl = qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        Cluster* best_a = nullptr;
        double best_a_ll = params.mu_log_alpha + next_ma.at(nullptr).ll;
        std::unordered_set<Cluster*>& matching_next_as = (
            xi[l+1] == -1 ? rs[l+1] : rs_by_emit[idx2d(l+1, xi[l+1], HP.K)]
        );
        for (Cluster *a : matching_next_as) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + get_msg_ll(next_ma, a);
            if (ll > best_a_ll) {
                best_a = a;
                best_a_ll = ll;
            }
        }
        double new_b_ll = -elogy + best_a_ll;
        mb[nullptr] = Msg{new_b_ll, best_a};

        // a messages.
        ma[nullptr] = Msg{new_a_ll + new_b_ll, nullptr};
        for (Cluster* a : matching_as) {
            Cluster* best_b = nullptr;
            double nFl = a->children.size();
            double best_b_ll = std::log(nFl) + params.mu_log_d[l] + mb[nullptr].ll;
            for (Cluster* b : a->children) {
                double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n) + get_msg_ll(mb, b);
                if (ll > best_b_ll) {
                    best_b = b;
                    best_b_ll = ll;
                }
            }
            best_b_ll -= std::log(static_cast<double>(a->n));
            ma[a] = Msg{best_b_ll, best_b};
        }
    }

    // Viterbi path.
    Cluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
        [](const auto& a, const auto& b) { return a.second.ll < b.second.ll; }
    )->first;
    int emission = xi[0] == -1 ? cluster_mode(0) : xi[0];
    Cluster* a_obj = (a == nullptr) ? create_cluster(std::vector<int>{}, true, 0, emission) : a;
    cluster_add(a_obj, i);

    Cluster* b = nullptr;
    Cluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = a_msgs[l].at(a).next;
        b_obj = (b == nullptr) ? create_cluster(std::vector<int>{}, false, l, -1) : b;
        cluster_add(b_obj, i);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = b_msgs[l].at(b).next;
        emission = xi[l+1] == -1 ? cluster_mode(l+1) : xi[l+1];
        a_obj = (a == nullptr) ? create_cluster(std::vector<int>{}, true, l+1, emission) : a;
        cluster_add(a_obj, i);
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
}

void Clusters::add_seqs(const std::vector<char>& x, const Params& params, HyperParams& HP) {
    int n = x.size() / HP.L;
    int old_N = HP.N;
    HP.N += n;
    r_assign.resize(r_assign.size() + x.size(), nullptr);
    q_assign.resize(q_assign.size() + n * (HP.L-1), nullptr);

    for (int i = 0; i < n; ++i) {
        viterbi_seq(x.begin() + i*HP.L, old_N + i, params);
    }
}

