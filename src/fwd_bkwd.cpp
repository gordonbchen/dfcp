#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "fwd_bkwd.hpp"
#include "seq_array.hpp"
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"
#include "max.hpp"
#include "math.hpp"
#include "util.hpp"


double get_msg_ll(const std::unordered_map<Cluster*, double>& msgs, Cluster* c, EmitMode emit_mode) {
    auto it = msgs.find(c);
    if (it == msgs.end()) {
        if (emit_mode != EmitMode::hard) {
            throw std::runtime_error("All msgs should be present unless emissions are hard.");
        }
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
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<std::unordered_map<Cluster*, double>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, double>> b_msgs(HP.L-1);
    for (int l = HP.L-1; l >= 0; --l) {
        int xil = get_xil(x, i, l, &obs_ls);
        double new_a_ll = get_cluster_emission_ll(nullptr, xil, l, clusters, params, HP);

        const std::unordered_set<Cluster*>& matching_as =
            (clusters.emit_mode != EmitMode::hard || xil == -1) ?
            clusters.rs[l] : clusters.rs_by_emit[idx2d(l, xil, HP.K)];
        if (l == HP.L-1) {
            a_msgs[l][nullptr] = new_a_ll;
            for (Cluster *a : matching_as) {
                a_msgs[l][a] = get_cluster_emission_ll(a, xil, l, clusters, params, HP);
            }
            continue;
        }

        // b messages.
        for (Cluster* b : clusters.qs[l]) {
            if (b->children.size() != 1) { throw std::runtime_error("b clusters can only have 1 child."); };
            Cluster* next_a = *b->children.begin();
            b_msgs[l][b] = get_msg_ll(a_msgs[l+1], next_a, clusters.emit_mode);
        }

        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        double new_b_ll = -elogy + params.mu_log_alpha + a_msgs[l+1].at(nullptr);
        int xil1 = get_xil(x, i, l+1, &obs_ls);
        const std::unordered_set<Cluster*>& matching_next_as =
            (clusters.emit_mode != EmitMode::hard || xil1 == -1)
            ? clusters.rs[l+1] : clusters.rs_by_emit[idx2d(l+1, xil1, HP.K)];
        for (Cluster *a : matching_next_as) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + get_msg_ll(a_msgs[l+1], a, clusters.emit_mode);
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
                    + get_msg_ll(b_msgs[l], b, clusters.emit_mode);
                a_ll = log_sum_exp(a_ll, ll);
            }
            a_msgs[l][a] = get_cluster_emission_ll(a, xil, l, clusters, params, HP) + a_ll;
        }
    }
    return a_msgs;
}

std::vector<std::unordered_map<Cluster*, double>> get_fwd_msgs(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<std::unordered_map<Cluster*, double>> a_msgs(HP.L);
    std::vector<std::unordered_map<Cluster*, double>> b_msgs(HP.L-1);
    for (int l = 0; l < HP.L; ++l) {
        int xil = get_xil(x, i, l, &obs_ls);
        double new_a_ll = get_cluster_emission_ll(nullptr, xil, l, clusters, params, HP);

        const std::unordered_set<Cluster*>& matching_as =
            (clusters.emit_mode != EmitMode::hard || xil == -1) ?
            clusters.rs[l] : clusters.rs_by_emit[idx2d(l, xil, HP.K)];
        if (l == 0) {
            a_msgs[l][nullptr] = new_a_ll + params.mu_log_alpha;
            for (Cluster *a : matching_as) {
                a_msgs[l][a] = get_cluster_emission_ll(a, xil, l, clusters, params, HP)
                    + std::log(static_cast<double>(a->n));
            }
            continue;
        }

        // b messages.
        for (Cluster* b : clusters.qs[l-1]) {
            if (b->parents.size() != 1) { throw std::runtime_error("b clusters can only have 1 parent."); };
            Cluster* prev_a = *b->parents.begin();
            double log_na = std::log(static_cast<double>(prev_a->n));
            double ab_ll = -log_na + delta_Elogx(params.mu_d[l-1], params.sigma2_d[l-1], -1.0, b->n);
            b_msgs[l-1][b] = get_msg_ll(a_msgs[l-1], prev_a, clusters.emit_mode) + ab_ll;
        }

        double new_b_ll = a_msgs[l-1].at(nullptr);
        int xilm1 = get_xil(x, i, l-1, &obs_ls);
        const std::unordered_set<Cluster*>& matching_prev_as =
            (clusters.emit_mode != EmitMode::hard || xilm1 == -1)
            ? clusters.rs[l-1] : clusters.rs_by_emit[idx2d(l-1, xilm1, HP.K)];
        for (Cluster *a : matching_prev_as) {
            double log_na = std::log(static_cast<double>(a->n));
            double nFl = a->children.size();
            double ab_ll = -log_na + std::log(nFl) + params.mu_log_d[l-1];
            new_b_ll = log_sum_exp(new_b_ll, get_msg_ll(a_msgs[l-1], a, clusters.emit_mode) + ab_ll);
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
            double a_ll = get_msg_ll(b_msgs[l-1], nullptr, clusters.emit_mode)
                -elogy + params.mu_log_d[l-1] + std::log(nCl);
            for (Cluster* b : a->parents) {
                a_ll = log_sum_exp(a_ll, get_msg_ll(b_msgs[l-1], b, clusters.emit_mode));
            }
            a_msgs[l][a] = get_cluster_emission_ll(a, xil, l, clusters, params, HP) + a_ll;
        }
    }
    return a_msgs;
}

void normalize_ll(std::vector<double>& ll, int L, int K) {
    for (int l = 0; l < L; ++l) {
        double sum_ll = ll[idx2d(l,0,K)];
        for (int k = 1; k < K; ++k) {
            sum_ll = log_sum_exp(sum_ll, ll[idx2d(l,k,K)]);
        }
        for (int k = 0; k < K; ++k) {
            ll[idx2d(l,k,K)] = std::exp(ll[idx2d(l,k,K)] - sum_ll);
        }
    }
}

std::vector<double> fwd_bkwd(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    auto a_msgs = get_bkwd_msgs(x, i, obs_ls, clusters, params, HP);
    auto fwd_a_msgs = get_fwd_msgs(x, i, obs_ls, clusters, params, HP);

    for (int l = 0; l < HP.L; ++l) {
        for (const auto& [a, ll] : fwd_a_msgs[l]) {
            a_msgs[l][a] = get_msg_ll(a_msgs[l], a, clusters.emit_mode) + ll;

            // Emission probs are double counted.
            int xil = get_xil(x, i, l, &obs_ls);
            a_msgs[l][a] -= get_cluster_emission_ll(a, xil, l, clusters, params, HP);
        }
    }

    int n_masked_ls = HP.L - obs_ls.size();
    std::vector<double> probs(n_masked_ls * HP.K, -std::numeric_limits<double>::infinity());
    int masked_l = 0;
    for (int l = 0; l < HP.L; ++l) {
        if (obs_ls.contains(l)) { continue; }
        for (const auto& [a, ll] : a_msgs[l]) {
            for (int k = 0; k < HP.K; ++k) {
                double emission_ll = get_cluster_emission_ll(a, k, l, clusters, params, HP);
                std::size_t prob_idx = idx2d(masked_l, k, HP.K);
                probs[prob_idx] = log_sum_exp(probs[prob_idx], ll + emission_ll);
            }
        }
        ++masked_l;
    }
    normalize_ll(probs, n_masked_ls, HP.K);
    return probs;
}
