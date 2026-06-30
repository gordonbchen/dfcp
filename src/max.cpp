#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "max.hpp"
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "soft_clusters.hpp"
#include "math.hpp"
#include "util.hpp"


double get_msg_ll(const std::unordered_map<Cluster*, Msg<Cluster>>& msgs, Cluster* c) {
    auto it = msgs.find(c);
    if (it == msgs.end()) {
        return -std::numeric_limits<double>::infinity();
    }
    return it->second.ll;
}

void viterbi_seq(Clusters& clusters, const std::vector<char>& x, int i, const HyperParams& HP, const Params& params) {
    std::vector<std::unordered_map<Cluster*, Msg<Cluster>>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, Msg<Cluster>>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        // Likelihood for new cluster.
        auto& ma = a_msgs[l];
        int emission = x[idx2d(i,l,HP.L)] == -1 ? clusters.cluster_mode(l) : x[idx2d(i,l,HP.L)];
        double new_a_ll = (
            delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, clusters.rs_by_emit[idx2d(l,emission,HP.K)].size())
            - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, clusters.rs[l].size())
        );

        std::unordered_set<Cluster*>& matching_as = (
            x[idx2d(i,l,HP.L)] == -1 ? clusters.rs[l] : clusters.rs_by_emit[idx2d(l, emission, HP.K)]
        );
        if (l == HP.L-1) {
            ma[nullptr] = Msg<Cluster>{new_a_ll, nullptr};
            for (Cluster *a : matching_as) {
                ma[a] = Msg<Cluster>{0.0, nullptr};
            }
            continue;
        }

        // b messages.
        auto& next_ma = a_msgs[l+1];
        auto& mb = b_msgs[l];
        for (Cluster* b : clusters.qs[l]) {
            if (b->children.size() != 1) { throw std::runtime_error("b clusters should only have 1 child."); };
            Cluster* next_a = *b->children.begin();
            mb[b] = Msg<Cluster>{get_msg_ll(next_ma, next_a), next_a};
         }

        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        Cluster* best_a = nullptr;
        double best_a_ll = params.mu_log_alpha + next_ma.at(nullptr).ll;
        std::unordered_set<Cluster*>& matching_next_as = (
            x[idx2d(i,l+1,HP.L)] == -1 ? clusters.rs[l+1] : clusters.rs_by_emit[idx2d(l+1, x[idx2d(i,l+1,HP.L)], HP.K)]
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
        mb[nullptr] = Msg<Cluster>{new_b_ll, best_a};

        // a messages.
        ma[nullptr] = Msg<Cluster>{new_a_ll + new_b_ll, nullptr};
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
            ma[a] = Msg<Cluster>{best_b_ll, best_b};
        }
    }

    // Viterbi path.
    Cluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
        [](const auto& a, const auto& b) { return a.second.ll < b.second.ll; }
    )->first;
    int emission = x[idx2d(i,0,HP.L)] == -1 ? clusters.cluster_mode(0) : x[idx2d(i,0,HP.L)];
    Cluster* a_obj = (a == nullptr) ? clusters.create_cluster(std::vector<int>{}, true, 0, emission) : a;
    clusters.cluster_add(a_obj, i);

    Cluster* b = nullptr;
    Cluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = a_msgs[l].at(a).next;
        b_obj = (b == nullptr) ? clusters.create_cluster(std::vector<int>{}, false, l, -1) : b;
        clusters.cluster_add(b_obj, i);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = b_msgs[l].at(b).next;
        emission = x[idx2d(i,l+1,HP.L)] == -1 ? clusters.cluster_mode(l+1) : x[idx2d(i,l+1,HP.L)];
        a_obj = (a == nullptr) ? clusters.create_cluster(std::vector<int>{}, true, l+1, emission) : a;
        clusters.cluster_add(a_obj, i);
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
}

void max_step(Clusters& clusters, const std::vector<char>& x, const HyperParams& HP, const Params& params) {
    for (int i = 0; i < HP.N; ++i) {
        for (int l = 0; l < HP.L; ++l) {
            clusters.cluster_remove(clusters.r_assign[idx2d(i, l, HP.L)], i);
            if (l == HP.L-1) {
                break;
            }
            clusters.cluster_remove(clusters.q_assign[idx2d(i, l, HP.L-1)], i);
        }
        viterbi_seq(clusters, x, i, HP, params);
    }
}


void viterbi_seq(
    SoftClusters& clusters, 
    const std::vector<char>& x, int i,
    const HyperParams& HP, const Params& params
) {
    std::vector<std::unordered_map<SoftCluster*, Msg<SoftCluster>>> a_msgs(HP.L);
    std::vector<std::unordered_map<SoftCluster*, Msg<SoftCluster>>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        // Likelihood for new cluster.
        auto& ma = a_msgs[l];
        double new_a_ll = 1.0 / HP.K;

        if (l == HP.L-1) {
            ma[nullptr] = Msg<SoftCluster>{new_a_ll, nullptr};
            for (SoftCluster *a : clusters.rs[l]) {
                int emission_count = (x[idx2d(i,l,HP.L)] == -1) ? a->mode().count : a->nk[x[idx2d(i,l,HP.L)]];
                double ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, emission_count);
                ll -= delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n);
                ma[a] = Msg<SoftCluster>{ll, nullptr};
            }
            continue;
        }

        // b messages.
        auto& next_ma = a_msgs[l+1];
        auto& mb = b_msgs[l];
        for (SoftCluster* b : clusters.qs[l]) {
            if (b->children.size() != 1) { throw std::runtime_error("b clusters should only have 1 child."); };
            SoftCluster* next_a = *b->children.begin();
            mb[b] = Msg<SoftCluster>{next_ma.at(next_a).ll, next_a};
         }

        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        SoftCluster* best_a = nullptr;
        double best_a_ll = params.mu_log_alpha + next_ma.at(nullptr).ll;
        for (SoftCluster *a : clusters.rs[l+1]) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + next_ma.at(a).ll;
            if (ll > best_a_ll) {
                best_a = a;
                best_a_ll = ll;
            }
        }
        double new_b_ll = -elogy + best_a_ll;
        mb[nullptr] = Msg<SoftCluster>{new_b_ll, best_a};

        // a messages.
        ma[nullptr] = Msg<SoftCluster>{new_a_ll + new_b_ll, nullptr};
        for (SoftCluster* a : clusters.rs[l]) {
            int emission_count = (x[idx2d(i,l,HP.L)] == -1) ? a->mode().count : a->nk[x[idx2d(i,l,HP.L)]];
            double emission_ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, emission_count);
            emission_ll -= delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n);

            SoftCluster* best_b = nullptr;
            double best_b_ll = std::log(a->children.size()) + params.mu_log_d[l] + mb[nullptr].ll;
            for (SoftCluster* b : a->children) {
                double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n) + mb.at(b).ll;
                if (ll > best_b_ll) {
                    best_b = b;
                    best_b_ll = ll;
                }
            }
            ma[a] = Msg<SoftCluster>{emission_ll - std::log(static_cast<double>(a->n)) + best_b_ll, best_b};
        }
    }

    // Viterbi path.
    SoftCluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
        [](const auto& a, const auto& b) { return a.second.ll < b.second.ll; }
    )->first;
    SoftCluster* a_obj = (a == nullptr) ? clusters.create_cluster(std::vector<int>{}, x, true, 0) : a;
    clusters.cluster_add(a_obj, i, x[idx2d(i,0,HP.L)]);

    SoftCluster* b = nullptr;
    SoftCluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = a_msgs[l].at(a).next;
        b_obj = (b == nullptr) ? clusters.create_cluster(std::vector<int>{}, x, false, l) : b;
        clusters.cluster_add(b_obj, i, -1);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = b_msgs[l].at(b).next;
        a_obj = (a == nullptr) ? clusters.create_cluster(std::vector<int>{}, x, true, l+1) : a;
        clusters.cluster_add(a_obj, i, x[idx2d(i,l,HP.L)]);
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
}

void max_step(SoftClusters& clusters, const std::vector<char>& x, const HyperParams& HP, const Params& params) {
    for (int i = 0; i < HP.N; ++i) {
        for (int l = 0; l < HP.L; ++l) {
            clusters.cluster_remove(clusters.r_assign[idx2d(i, l, HP.L)], i, x[idx2d(i, l, HP.L)]);
            if (l == HP.L-1) {
                break;
            }
            clusters.cluster_remove(clusters.q_assign[idx2d(i, l, HP.L-1)], i, -1);
        }
        viterbi_seq(clusters, x, i, HP, params);
    }
}

