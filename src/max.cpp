#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <boost/math/special_functions/digamma.hpp>
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

int get_new_cluster_emission(
    int8_t xil, int l,
    double Elog_match, double Elog_mismatch,
    const Clusters& clusters, const HyperParams& HP, const Params& params
) {
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
            + (xil == k ? Elog_match : Elog_mismatch);
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
    double Elog_match, Elog_mismatch, mu_eps, sigma2_eps;
    Elog_match = Elog_mismatch = mu_eps = sigma2_eps = std::numeric_limits<double>::infinity();
    if (clusters.noisy) {
        double digamma_eps_sum = boost::math::digamma(params.alpha_eps + params.beta_eps);
        Elog_match = boost::math::digamma(params.beta_eps) - digamma_eps_sum;
        Elog_mismatch = boost::math::digamma(params.alpha_eps) - digamma_eps_sum - std::log(HP.K-1.0);
        mu_eps = params.alpha_eps / (params.alpha_eps + params.beta_eps);
        sigma2_eps = (params.alpha_eps * params.beta_eps)
            / (std::pow(params.alpha_eps + params.beta_eps, 2.0) * (params.alpha_eps + params.beta_eps + 1.0));
    }

    std::vector<std::unordered_map<Cluster*, Msg>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, Msg>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        // Likelihood for new cluster.
        double new_a_ll = 0.0;
        if (xi[l] != -1) {
            if (clusters.soft) {
                new_a_ll = -std::log(static_cast<double>(HP.K));
            }
            else {
                int nkl = clusters.rs_by_emit[idx2d(l,xi[l],HP.K)].size();
                new_a_ll = -delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, clusters.rs[l].size(), params.mu_log_gamma[l]);
                if (clusters.noisy) {
                    double c = (static_cast<double>(clusters.rs[l].size()) - HP.K*nkl) / (HP.K-1.0);
                    double mu_y = params.mu_gamma[l] + nkl + c*mu_eps;
                    double sigma2_y = params.sigma2_gamma[l] + c*c*sigma2_eps;
                    new_a_ll += delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);
                }
                else {
                    new_a_ll += delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nkl, params.mu_log_gamma[l]);
                }
            }
        }

        std::unordered_set<Cluster*>& matching_as = clusters.noisy || clusters.soft || (xi[l] == -1) ?
            clusters.rs[l] : clusters.rs_by_emit[idx2d(l, xi[l], HP.K)];
        if (l == HP.L-1) {
            a_msgs[l][nullptr] = Msg{new_a_ll, nullptr};
            for (Cluster *a : matching_as) {
                double ll = 0.0;
                if (clusters.noisy && xi[l] != -1) {
                    ll = xi[l] == a->emission ? Elog_match : Elog_mismatch;
                }
                else if (clusters.soft && xi[l] != -1) {
                    ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, a->nk[xi[l]], params.mu_log_gamma[l])
                        - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n_obs, params.mu_log_gamma[l]);
                }
                a_msgs[l][a] = Msg{ll, nullptr};
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
        std::unordered_set<Cluster*>& matching_next_as = clusters.noisy || clusters.soft || (xi[l+1] == -1) ?
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
            double emission_ll = 0.0;
            if (clusters.noisy && (xi[l] != -1)) {
                emission_ll = xi[l] == a->emission ? Elog_match : Elog_mismatch;
            }
            if (clusters.soft && (xi[l] != -1)) {
                emission_ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, a->nk[xi[l]], params.mu_log_gamma[l])
                    - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n_obs, params.mu_log_gamma[l]);
            }
            a_msgs[l][a] = Msg{emission_ll - std::log(static_cast<double>(a->n)) + best_b_ll, best_b};
        }
    }

    // Viterbi path.
    Cluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
        [](const auto& a, const auto& b) { return a.second.ll < b.second.ll; }
    )->first;
    Cluster* a_obj = a;
    if (a == nullptr) {
        int emission = get_new_cluster_emission(xi[0], 0, Elog_match, Elog_mismatch, clusters, HP, params);
        a_obj = clusters.create_empty_cluster(true, 0, emission);
    }
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
        a_obj = a;
        if (a == nullptr) {
            int emission = get_new_cluster_emission(xi[l+1], l+1, Elog_match, Elog_mismatch, clusters, HP, params);
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

    double digamma_eps_sum = boost::math::digamma(params.alpha_eps + params.beta_eps);
    double Elog_match = boost::math::digamma(params.beta_eps) - digamma_eps_sum;
    double Elog_mismatch = boost::math::digamma(params.alpha_eps) - digamma_eps_sum - std::log(HP.K-1.0);

    for (int l = 0; l < HP.L; ++l) {
        for (Cluster* a : clusters.rs[l]) {
            int best_k = 0;
            double best_ll = -std::numeric_limits<double>::infinity();
            for (int k = 0; k < HP.K; ++k) {
                int nkl = clusters.rs_by_emit[idx2d(l,k,HP.K)].size() - (a->emission == k);
                double ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nkl, params.mu_log_gamma[l])
                    + a->nk[k]*Elog_match + (a->n_obs - a->nk[k])*Elog_mismatch;
                if (ll > best_ll) {
                    best_ll = ll;
                    best_k = k;
                }
            }
            clusters.set_emission(a, best_k);
        }
    }
}

