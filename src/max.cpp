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

double get_msg_ll(const std::unordered_map<Cluster*, Msg>& msgs, Cluster* c, bool noisy, bool soft) {
    auto it = msgs.find(c);
    if (it == msgs.end()) {
        if (noisy || soft) { throw std::runtime_error("All msgs should be present if noisy or soft."); }
        return -std::numeric_limits<double>::infinity();
    }
    return it->second.ll;
}

double get_cluster_emission_ll(
    Cluster* a, int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    // TODO: noisy, soft, hard enum.
    if (xil == -1) {
        return 0.0;
    }

    if (a == nullptr) {
        if (clusters.soft) {
            return -std::log(static_cast<double>(HP.K));
        }

        int nkl = clusters.rs_by_emit[idx2d(l,xil,HP.K)].size();
        double ll = -delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, clusters.rs[l].size(), params.mu_log_gamma[l]);
        if (clusters.noisy) {
            double c = (static_cast<double>(clusters.rs[l].size()) - HP.K*nkl) / (HP.K-1.0);
            double mu_y = params.mu_gamma[l] + nkl + c*params.mu_eps;
            double sigma2_y = params.sigma2_gamma[l] + c*c*params.sigma2_eps;
            ll += delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);
            return ll;
        }

        ll += delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nkl, params.mu_log_gamma[l]);
        return ll;
    }

    if (clusters.soft) {
        return delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, a->nk[xil], params.mu_log_gamma[l])
            - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n_obs, params.mu_log_gamma[l]);
    }
    if (clusters.noisy) {
        return xil == a->emission ? params.Eeps_log_match : params.Eeps_log_mismatch;
    }
    return xil == a->emission ? 0.0 : -std::numeric_limits<double>::infinity();
}

std::vector<Cluster*> get_viterbi_clusters(
    const Clusters& clusters, std::vector<int8_t>::const_iterator xi,
    const HyperParams& HP, const Params& params
) {
    std::vector<std::unordered_map<Cluster*, Msg>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, Msg>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        double new_a_ll = get_cluster_emission_ll(nullptr, xi[l], l, clusters, params, HP);

        const std::unordered_set<Cluster*>& matching_as = clusters.noisy || clusters.soft || (xi[l] == -1) ?
            clusters.rs[l] : clusters.rs_by_emit[idx2d(l, xi[l], HP.K)];
        if (l == HP.L-1) {
            a_msgs[l][nullptr] = Msg{new_a_ll, nullptr};
            for (Cluster *a : matching_as) {
                a_msgs[l][a] = Msg{get_cluster_emission_ll(a, xi[l], l, clusters, params, HP), nullptr};
            }
            continue;
        }

        // b messages.
        for (Cluster* b : clusters.qs[l]) {
            if (b->children.size() != 1) { throw std::runtime_error("b clusters should only have 1 child."); };
            Cluster* next_a = *b->children.begin();
            b_msgs[l][b] = Msg{get_msg_ll(a_msgs[l+1], next_a, clusters.noisy, clusters.soft), next_a};
        }

        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        Cluster* best_a = nullptr;
        double best_a_ll = params.mu_log_alpha + a_msgs[l+1].at(nullptr).ll;
        const std::unordered_set<Cluster*>& matching_next_as = clusters.noisy || clusters.soft || (xi[l+1] == -1) ?
            clusters.rs[l+1] : clusters.rs_by_emit[idx2d(l+1, xi[l+1], HP.K)];
        for (Cluster *a : matching_next_as) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + get_msg_ll(a_msgs[l+1], a, clusters.noisy, clusters.soft);
            if (ll > best_a_ll) {
                best_a = a;
                best_a_ll = ll;
            }
        }
        double new_b_ll = -elogy + best_a_ll;
        b_msgs[l][nullptr] = Msg{new_b_ll, best_a};

        // a messages.
        a_msgs[l][nullptr] = Msg{new_a_ll + new_b_ll, nullptr};
        for (Cluster* a : matching_as) {
            Cluster* best_b = nullptr;
            double nFl = a->children.size();
            double best_b_ll = std::log(nFl) + params.mu_log_d[l] + b_msgs[l][nullptr].ll;
            for (Cluster* b : a->children) {
                double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n)
                    + get_msg_ll(b_msgs[l], b, clusters.noisy, clusters.soft);
                if (ll > best_b_ll) {
                    best_b = b;
                    best_b_ll = ll;
                }
            }
            double emission_ll = get_cluster_emission_ll(a, xi[l], l, clusters, params, HP);
            a_msgs[l][a] = Msg{emission_ll - std::log(static_cast<double>(a->n)) + best_b_ll, best_b};
        }
    }

    // Initial CRP probs.
    for (auto& [a, m] : a_msgs[0]) {
        m.ll += (a == nullptr) ? params.mu_log_alpha : std::log(static_cast<double>(a->n));
    }

    // Viterbi path.
    std::vector<Cluster*> best_clusters;
    best_clusters.reserve(HP.L + HP.L-1);

    Cluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
        [](const auto& a, const auto& b) { return a.second.ll < b.second.ll; }
    )->first;
    best_clusters.emplace_back(a);

    Cluster* b;
    for (int l = 0; l < HP.L-1; ++l) {
        b = a_msgs[l].at(a).next;
        best_clusters.emplace_back(b);

        a = b_msgs[l].at(b).next;
        best_clusters.emplace_back(a);
    }
    return best_clusters;
}


int get_new_cluster_emission(int8_t xil, int l, const Clusters& clusters, const HyperParams& HP, const Params& params) {
    if (clusters.soft) {
        return -1;
    }
    if (xil == -1) {
        return clusters.cluster_mode(l);
    }
    if (!clusters.noisy) {
        return xil;
    }

    int best_k = 0;
    double best_ll = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < HP.K; ++k) {
        int nkl = clusters.rs_by_emit[idx2d(l,k,HP.K)].size();
        double ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nkl, params.mu_log_gamma[l])
            + (xil == k ? params.Eeps_log_match : params.Eeps_log_mismatch);
        if (ll > best_ll) {
            best_ll = ll;
            best_k = k;
        }
    }
    return best_k;
}

void viterbi_seq(
    Clusters& clusters, std::vector<int8_t>::const_iterator xi, int i,
    const HyperParams& HP, const Params& params
) {
    std::vector<Cluster*> best_clusters{get_viterbi_clusters(clusters, xi, HP, params)};

    Cluster* a = best_clusters[0];
    Cluster* a_obj = a;
    if (a == nullptr) {
        int emission = get_new_cluster_emission(xi[0], 0, clusters, HP, params);
        a_obj = clusters.create_empty_cluster(true, 0, emission);
    }
    clusters.cluster_add(a_obj, i, xi[0]);

    Cluster* b = nullptr;
    Cluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = best_clusters[2*l + 1];
        b_obj = (b == nullptr) ? clusters.create_empty_cluster(false, l, -1) : b;
        clusters.cluster_add(b_obj, i, -1);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = best_clusters[2 * (l+1)];
        a_obj = a;
        if (a == nullptr) {
            int emission = get_new_cluster_emission(xi[l+1], l+1, clusters, HP, params);
            a_obj = clusters.create_empty_cluster(true, l+1, emission);
        }
        clusters.cluster_add(a_obj, i, xi[l+1]);
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
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


void max_cluster_emissions(Clusters& clusters, const HyperParams& HP, const Params& params) {
    if (!clusters.noisy) { throw std::runtime_error("Only need to maximize cluster emissions if noisy."); }
    for (int l = 0; l < HP.L; ++l) {
        for (Cluster* a : clusters.rs[l]) {
            int best_k = 0;
            double best_ll = -std::numeric_limits<double>::infinity();
            for (int k = 0; k < HP.K; ++k) {
                int nkl = clusters.rs_by_emit[idx2d(l,k,HP.K)].size() - (a->emission == k);
                double ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nkl, params.mu_log_gamma[l])
                    + a->nk[k]*params.Eeps_log_match + (a->n_obs - a->nk[k])*params.Eeps_log_mismatch;
                if (ll > best_ll) {
                    best_ll = ll;
                    best_k = k;
                }
            }
            clusters.set_emission(a, best_k);
        }
    }
}

