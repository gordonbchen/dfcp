#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "max.hpp"
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "math.hpp"
#include "util.hpp"


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

void hard_viterbi_seq(
    Clusters& clusters, std::vector<int8_t>::const_iterator xi, int i,
    const HyperParams& HP, const Params& params
) {
    std::vector<std::unordered_map<Cluster*, Msg>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, Msg>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        // Likelihood for new cluster.
        auto& ma = a_msgs[l];
        int emission = xi[l] == -1 ? clusters.cluster_mode(l) : xi[l];
        int nkl = clusters.rs_by_emit[idx2d(l,emission,HP.K)].size();
        double new_a_ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nkl);
        new_a_ll -= delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, clusters.rs[l].size());

        std::unordered_set<Cluster*>& matching_as = (
            xi[l] == -1 ? clusters.rs[l] : clusters.rs_by_emit[idx2d(l, emission, HP.K)]
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
        for (Cluster* b : clusters.qs[l]) {
            if (b->children.size() != 1) {
                throw std::runtime_error("b clusters should only have 1 child.");
            };
            Cluster* next_a = *b->children.begin();
            mb[b] = Msg{get_msg_ll(next_ma, next_a), next_a};
         }

        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        Cluster* best_a = nullptr;
        double best_a_ll = params.mu_log_alpha + next_ma.at(nullptr).ll;
        std::unordered_set<Cluster*>& matching_next_as = (
            xi[l+1] == -1 ? clusters.rs[l+1] : clusters.rs_by_emit[idx2d(l+1, xi[l+1], HP.K)]
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
    int emission = xi[0] == -1 ? clusters.cluster_mode(0) : xi[0];
    Cluster* a_obj = (a == nullptr) ? clusters.create_empty_cluster(true, 0, emission) : a;
    clusters.cluster_add(a_obj, i, emission);

    Cluster* b = nullptr;
    Cluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = a_msgs[l].at(a).next;
        b_obj = (b == nullptr) ? clusters.create_empty_cluster(false, l, -1) : b;
        clusters.cluster_add(b_obj, i, -1);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = b_msgs[l].at(b).next;
        emission = xi[l+1] == -1 ? clusters.cluster_mode(l+1) : xi[l+1];
        a_obj = (a == nullptr) ? clusters.create_empty_cluster(true, l+1, emission) : a;
        clusters.cluster_add(a_obj, i, emission);
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
}

void soft_viterbi_seq(
    Clusters& clusters, std::vector<int8_t>::const_iterator xi, int i,
    const HyperParams& HP, const Params& params
) {
    std::vector<std::unordered_map<Cluster*, Msg>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, Msg>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        // Likelihood for new cluster.
        auto& ma = a_msgs[l];
        double new_a_ll = xi[l] == -1 ? 0.0 : -std::log(static_cast<double>(HP.K));

        if (l == HP.L-1) {
            ma[nullptr] = Msg{new_a_ll, nullptr};
            for (Cluster *a : clusters.rs[l]) {
                double ll = 0.0;
                if (xi[l] != -1) {
                    int emission_count = a->nk[xi[l]];
                    ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, emission_count);
                    ll -= delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n_obs);
                }
                ma[a] = Msg{ll, nullptr};
            }
            continue;
        }

        // b messages.
        auto& next_ma = a_msgs[l+1];
        auto& mb = b_msgs[l];
        for (Cluster* b : clusters.qs[l]) {
            if (b->children.size() != 1) {
                throw std::runtime_error("b clusters should only have 1 child.");
            };
            Cluster* next_a = *b->children.begin();
            mb[b] = Msg{next_ma.at(next_a).ll, next_a};
         }

        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        Cluster* best_a = nullptr;
        double best_a_ll = params.mu_log_alpha + next_ma.at(nullptr).ll;
        for (Cluster *a : clusters.rs[l+1]) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + next_ma.at(a).ll;
            if (ll > best_a_ll) {
                best_a = a;
                best_a_ll = ll;
            }
        }
        double new_b_ll = -elogy + best_a_ll;
        mb[nullptr] = Msg{new_b_ll, best_a};

        // a messages.
        ma[nullptr] = Msg{new_a_ll + new_b_ll, nullptr};
        for (Cluster* a : clusters.rs[l]) {
            double emission_ll = 0.0;
            if (xi[l] != -1) {
                int emission_count = a->nk[xi[l]];
                emission_ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, emission_count);
                emission_ll -= delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n_obs);
            }

            Cluster* best_b = nullptr;
            double nFl = a->children.size();
            double best_b_ll = std::log(nFl) + params.mu_log_d[l] + mb[nullptr].ll;
            for (Cluster* b : a->children) {
                double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n) + mb.at(b).ll;
                if (ll > best_b_ll) {
                    best_b = b;
                    best_b_ll = ll;
                }
            }
            ma[a] = Msg{emission_ll - std::log(static_cast<double>(a->n)) + best_b_ll, best_b};
        }
    }

    // Viterbi path.
    Cluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
        [](const auto& a, const auto& b) { return a.second.ll < b.second.ll; }
    )->first;
    Cluster* a_obj = (a == nullptr) ? clusters.create_empty_cluster(true, 0, -1) : a;
    clusters.cluster_add(a_obj, i, xi[0]);

    Cluster* b = nullptr;
    Cluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = a_msgs[l].at(a).next;
        b_obj = (b == nullptr) ? clusters.create_empty_cluster(false, l, -1) : b;
        clusters.cluster_add(b_obj, i, -1);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = b_msgs[l].at(b).next;
        a_obj = (a == nullptr) ? clusters.create_empty_cluster(true, l+1, -1) : a;
        clusters.cluster_add(a_obj, i, xi[l+1]);
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
}

void viterbi_seq(
    Clusters& clusters, std::vector<int8_t>::const_iterator xi, int i,
    const HyperParams& HP, const Params& params
) {
    if (clusters.soft) {
        return soft_viterbi_seq(clusters, xi, i, HP, params);
    }
    return hard_viterbi_seq(clusters, xi, i, HP, params);
}

void max_step(
    Clusters& clusters, const std::vector<int8_t>& x, const HyperParams& HP, const Params& params
) {
    for (int i = 0; i < HP.N; ++i) {
        for (int l = 0; l < HP.L; ++l) {
            clusters.cluster_remove(clusters.r_assign[idx2d(i, l, HP.L)], i, x[idx2d(i, l, HP.L)]);
            if (l == HP.L-1) {
                break;
            }
            clusters.cluster_remove(clusters.q_assign[idx2d(i, l, HP.L-1)], i, -1);
        }
        viterbi_seq(clusters, x.begin() + i*HP.L, i, HP, params);
    }
}

void add_seqs(
    Clusters& clusters, std::vector<int8_t>::const_iterator new_x, int n_new,
    HyperParams& HP, const Params& params
) {
    int old_N = HP.N;
    HP.N += n_new;
    clusters.r_assign.resize(HP.N * HP.L, nullptr);
    clusters.q_assign.resize(HP.N * (HP.L-1), nullptr);

    for (int i = 0; i < n_new; ++i) {
        viterbi_seq(clusters, new_x + i*HP.L, old_N + i, HP, params);
    }
}

