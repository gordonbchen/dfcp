#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"
#include "max.hpp"
#include "math.hpp"
#include "util.hpp"


double get_msg_ll(const std::unordered_map<Cluster*, double>& msgs, Cluster* c, bool noisy, bool soft) {
    auto it = msgs.find(c);
    if (it == msgs.end()) {
        if (noisy || soft) { throw std::runtime_error("All msgs should be present if noisy or soft."); }
        return -std::numeric_limits<double>::infinity();
    }
    return it->second;
}

double log_sum_exp(double x1, double x2) {
    if (x1 == -std::numeric_limits<double>::infinity()) {
        return x2;
    }
    if (x2 == -std::numeric_limits<double>::infinity()) {
        return x1;
    }
    double xmax = std::max(x1, x2);
    return xmax + std::log(std::exp(x1 - xmax) + std::exp(x2 - xmax));
}

std::vector<std::unordered_map<Cluster*, double>> get_bkwd_msgs(
    std::vector<int8_t>::const_iterator xi, const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<std::unordered_map<Cluster*, double>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, double>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        double new_a_ll = get_new_cluster_ll(xi[l], l, clusters, params, HP);

        const std::unordered_set<Cluster*>& matching_as = clusters.noisy || clusters.soft || (xi[l] == -1) ?
            clusters.rs[l] : clusters.rs_by_emit[idx2d(l, xi[l], HP.K)];
        if (l == HP.L-1) {
            a_msgs[l][nullptr] = new_a_ll;
            for (Cluster *a : matching_as) {
                a_msgs[l][a] = get_cluster_emission_ll(a, xi[l], l, clusters, params, HP);
            }
            continue;
        }

        // b messages.
        for (Cluster* b : clusters.qs[l]) {
            if (b->children.size() != 1) { throw std::runtime_error("b clusters should only have 1 child."); };
            Cluster* next_a = *b->children.begin();
            b_msgs[l][b] = get_msg_ll(a_msgs[l+1], next_a, clusters.noisy, clusters.soft);
        }

        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        double new_b_ll = -elogy + params.mu_log_alpha + a_msgs[l+1].at(nullptr);
        const std::unordered_set<Cluster*>& matching_next_as = clusters.noisy || clusters.soft || (xi[l+1] == -1) ?
            clusters.rs[l+1] : clusters.rs_by_emit[idx2d(l+1, xi[l+1], HP.K)];
        for (Cluster *a : matching_next_as) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + get_msg_ll(a_msgs[l+1], a, clusters.noisy, clusters.soft);
            new_b_ll = log_sum_exp(new_b_ll, -elogy + ll);
        }
        b_msgs[l][nullptr] = new_b_ll;

        // a messages.
        a_msgs[l][nullptr] = new_a_ll + new_b_ll;
        for (Cluster* a : matching_as) {
            double log_na = std::log(static_cast<double>(a->n));
            double nFl = a->children.size();
            double a_ll = -log_na + std::log(nFl) + params.mu_log_d[l] + b_msgs[l][nullptr];
            for (Cluster* b : a->children) {
                double ll = -log_na + delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n)
                    + get_msg_ll(b_msgs[l], b, clusters.noisy, clusters.soft);
                a_ll = log_sum_exp(a_ll, ll);
            }
            a_msgs[l][a] = get_cluster_emission_ll(a, xi[l], l, clusters, params, HP) + a_ll;
        }
    }
    return a_msgs;
}

std::vector<std::unordered_map<Cluster*, double>> get_fwd_msgs(
    std::vector<int8_t>::const_iterator xi, const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<std::unordered_map<Cluster*, double>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, double>> b_msgs(HP.L-1);
    for (int l = 0; l < HP.L; ++l) {
        double new_a_ll = get_new_cluster_ll(xi[l], l, clusters, params, HP);

        const std::unordered_set<Cluster*>& matching_as = clusters.noisy || clusters.soft || (xi[l] == -1) ?
            clusters.rs[l] : clusters.rs_by_emit[idx2d(l, xi[l], HP.K)];
        if (l == 0) {
            a_msgs[l][nullptr] = new_a_ll + params.mu_log_alpha;
            for (Cluster *a : matching_as) {
                a_msgs[l][a] = get_cluster_emission_ll(a, xi[l], l, clusters, params, HP)
                    + std::log(static_cast<double>(a->n));
            }
            continue;
        }

        // b messages.
        for (Cluster* b : clusters.qs[l-1]) {
            if (b->parents.size() != 1) { throw std::runtime_error("b clusters should only have 1 parent."); };
            Cluster* prev_a = *b->parents.begin();
            double log_na = std::log(static_cast<double>(prev_a->n));
            double ab_ll = -log_na + delta_Elogx(params.mu_d[l-1], params.sigma2_d[l-1], -1.0, b->n);
            b_msgs[l-1][b] = get_msg_ll(a_msgs[l-1], prev_a, clusters.noisy, clusters.soft) + ab_ll;
        }

        double new_b_ll = a_msgs[l-1].at(nullptr);
        const std::unordered_set<Cluster*>& matching_prev_as = clusters.noisy || clusters.soft || (xi[l-1] == -1) ?
            clusters.rs[l-1] : clusters.rs_by_emit[idx2d(l-1, xi[l-1], HP.K)];
        for (Cluster *a : matching_prev_as) {
            double log_na = std::log(static_cast<double>(a->n));
            double nFl = a->children.size();
            double ab_ll = -log_na + std::log(nFl) + params.mu_log_d[l-1];
            new_b_ll = log_sum_exp(new_b_ll, get_msg_ll(a_msgs[l-1], a, clusters.noisy, clusters.soft) + ab_ll);
        }
        b_msgs[l-1][nullptr] = new_b_ll;

        // a messages.
        int nQl = clusters.qs[l-1].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l-1];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l-1];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);
        a_msgs[l][nullptr] = new_b_ll - elogy + params.mu_log_alpha + new_a_ll;

        for (Cluster* a : matching_as) {
            double nCl = a->parents.size();
            double a_ll = get_msg_ll(b_msgs[l-1], nullptr, clusters.noisy, clusters.soft)
                -elogy + params.mu_log_d[l-1] + std::log(nCl);
            for (Cluster* b : a->parents) {
                a_ll = log_sum_exp(a_ll, get_msg_ll(b_msgs[l-1], b, clusters.noisy, clusters.soft));
            }
            a_msgs[l][a] = get_cluster_emission_ll(a, xi[l], l, clusters, params, HP) + a_ll;
        }
    }
    return a_msgs;
}

std::vector<double> forward_backward(
    std::vector<int8_t>::const_iterator xi, const std::vector<int>& prob_idxs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<std::unordered_map<Cluster*, double>> a_msgs = get_bkwd_msgs(xi, clusters, params, HP);
    std::vector<std::unordered_map<Cluster*, double>> fwd_a_msgs = get_fwd_msgs(xi, clusters, params, HP);

    for (int l = 0; l < HP.L; ++l) {
        for (const auto& [a, ll] : fwd_a_msgs[l]) {
            a_msgs[l][a] = get_msg_ll(a_msgs[l], a, clusters.noisy, clusters.soft) + ll;
            // Emission probs are double counted.
            a_msgs[l][a] -= (a == nullptr) ? get_new_cluster_ll(xi[l], l, clusters, params, HP)
                : get_cluster_emission_ll(a, xi[l], l, clusters, params, HP);
        }
    }

    std::vector<double> probs(prob_idxs.size() * HP.K, -std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < prob_idxs.size(); ++i) {
        int l = prob_idxs[i];
        for (const auto& [a, ll] : a_msgs[l]) {
            for (int k = 0; k < HP.K; ++k) {
                double emission_ll = (a == nullptr) ? get_new_cluster_ll(k, l, clusters, params, HP)
                    : get_cluster_emission_ll(a, k, l, clusters, params, HP);
                probs[idx2d(i,k,HP.K)] = log_sum_exp(probs[idx2d(i,k,HP.K)], ll + emission_ll);
            }
        }

        // Normalize.
        double sum_ll = probs[idx2d(i,0,HP.K)];
        for (int k = 1; k < HP.K; ++k) {
            sum_ll = log_sum_exp(sum_ll, probs[idx2d(i,k,HP.K)]);
        }
        for (int k = 0; k < HP.K; ++k) {
            probs[idx2d(i,k,HP.K)] = std::exp(probs[idx2d(i,k,HP.K)] - sum_ll);
        }
    }
    return probs;
}

