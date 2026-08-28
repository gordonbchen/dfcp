#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>
#include "fwd_bkwd.hpp"
#include "seq_array.hpp"
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"
#include "max.hpp"
#include "math.hpp"
#include "util.hpp"


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

struct FwdBkwdMsgs {
    std::vector<double> a_msgs;
    std::vector<double> new_a_msgs;

    FwdBkwdMsgs(uint32_t n_cluster_ids, int L) : a_msgs(n_cluster_ids), new_a_msgs(L) {}
};

FwdBkwdMsgs get_bkwd_msgs(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    FwdBkwdMsgs msgs(clusters.next_cluster_id, HP.L);
    std::vector<double>& a_msgs = msgs.a_msgs;
    std::vector<double>& new_a_msgs = msgs.new_a_msgs;
    for (int l = HP.L-1; l >= 0; --l) {
        int xil = get_xil(x, i, l, &obs_ls);
        double new_a_ll = get_cluster_emission_ll(nullptr, xil, l, clusters, params, HP);

        const std::vector<Cluster*>& matching_as = clusters.get_matching_as(l, xil);
        if (l == HP.L-1) {
            new_a_msgs[l] = new_a_ll;
            for (Cluster *a : matching_as) {
                a_msgs[a->id] = get_cluster_emission_ll(a, xil, l, clusters, params, HP);
            }
            continue;
        }

        // New b message.
        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        double new_b_ll = -elogy + params.mu_log_alpha + new_a_msgs[l+1];
        int xil1 = get_xil(x, i, l+1, &obs_ls);
        for (Cluster *a : clusters.get_matching_as(l+1, xil1)) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + a_msgs[a->id];
            new_b_ll = log_sum_exp(new_b_ll, -elogy + ll);
        }

        // a messages.
        new_a_msgs[l] = new_a_ll + new_b_ll;
        for (Cluster* a : matching_as) {
            double log_na = std::log(static_cast<double>(a->n));
            double nFl = a->children.size();
            double a_ll = -log_na + std::log(nFl) + params.mu_log_d[l] + new_b_ll;
            for (Cluster* b : a->children) {
                Cluster* next_a = b->children[0];
                if (clusters.emit_mode == EmitMode::hard && xil1 != -1 && next_a->emission != xil1) {
                    continue;
                }
                double ll = -log_na + delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n)
                    + a_msgs[next_a->id];
                a_ll = log_sum_exp(a_ll, ll);
            }
            a_msgs[a->id] = get_cluster_emission_ll(a, xil, l, clusters, params, HP) + a_ll;
        }
    }
    return msgs;
}

FwdBkwdMsgs get_fwd_msgs(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    FwdBkwdMsgs msgs(clusters.next_cluster_id, HP.L);
    std::vector<double>& a_msgs = msgs.a_msgs;
    std::vector<double>& new_a_msgs = msgs.new_a_msgs;
    for (int l = 0; l < HP.L; ++l) {
        int xil = get_xil(x, i, l, &obs_ls);
        double new_a_ll = get_cluster_emission_ll(nullptr, xil, l, clusters, params, HP);

        const std::vector<Cluster*>& matching_as = clusters.get_matching_as(l, xil);
        if (l == 0) {
            new_a_msgs[l] = new_a_ll + params.mu_log_alpha;
            for (Cluster *a : matching_as) {
                a_msgs[a->id] = get_cluster_emission_ll(a, xil, l, clusters, params, HP)
                    + std::log(static_cast<double>(a->n));
            }
            continue;
        }

        // New b message.
        double new_b_ll = new_a_msgs[l-1];
        int xilm1 = get_xil(x, i, l-1, &obs_ls);
        for (Cluster *a : clusters.get_matching_as(l-1, xilm1)) {
            double log_na = std::log(static_cast<double>(a->n));
            double nFl = a->children.size();
            double ab_ll = -log_na + std::log(nFl) + params.mu_log_d[l-1];
            new_b_ll = log_sum_exp(new_b_ll, a_msgs[a->id] + ab_ll);
        }

        // a messages.
        int nQl = clusters.qs[l-1].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l-1];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l-1];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);
        new_a_msgs[l] = new_b_ll - elogy + params.mu_log_alpha + new_a_ll;

        for (Cluster* a : matching_as) {
            double nCl = a->parents.size();
            double a_ll = new_b_ll - elogy + params.mu_log_d[l-1] + std::log(nCl);
            for (Cluster* b : a->parents) {
                Cluster* prev_a = b->parents[0];
                if (clusters.emit_mode == EmitMode::hard && xilm1 != -1 && prev_a->emission != xilm1) {
                    continue;
                }
                double log_na = std::log(static_cast<double>(prev_a->n));
                double ab_ll = -log_na + delta_Elogx(params.mu_d[l-1], params.sigma2_d[l-1], -1.0, b->n);
                a_ll = log_sum_exp(a_ll, a_msgs[prev_a->id] + ab_ll);
            }
            a_msgs[a->id] = get_cluster_emission_ll(a, xil, l, clusters, params, HP) + a_ll;
        }
    }
    return msgs;
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
    FwdBkwdMsgs msgs = get_bkwd_msgs(x, i, obs_ls, clusters, params, HP);
    FwdBkwdMsgs fwd_msgs = get_fwd_msgs(x, i, obs_ls, clusters, params, HP);

    for (int l = 0; l < HP.L; ++l) {
        int xil = get_xil(x, i, l, &obs_ls);
        msgs.new_a_msgs[l] += fwd_msgs.new_a_msgs[l]
            - get_cluster_emission_ll(nullptr, xil, l, clusters, params, HP);
        for (Cluster* a : clusters.get_matching_as(l, xil)) {
            msgs.a_msgs[a->id] += fwd_msgs.a_msgs[a->id]
                - get_cluster_emission_ll(a, xil, l, clusters, params, HP);
        }
    }

    int n_masked_ls = HP.L - obs_ls.size();
    std::vector<double> probs(n_masked_ls * HP.K, -std::numeric_limits<double>::infinity());
    int masked_l = 0;
    for (int l = 0; l < HP.L; ++l) {
        if (obs_ls.contains(l)) { continue; }
        for (int k = 0; k < HP.K; ++k) {
            std::size_t prob_idx = idx2d(masked_l, k, HP.K);
            double emission_ll = get_cluster_emission_ll(nullptr, k, l, clusters, params, HP);
            probs[prob_idx] = log_sum_exp(probs[prob_idx], msgs.new_a_msgs[l] + emission_ll);
        }
        for (Cluster* a : clusters.rs[l]) {
            for (int k = 0; k < HP.K; ++k) {
                std::size_t prob_idx = idx2d(masked_l, k, HP.K);
                double emission_ll = get_cluster_emission_ll(a, k, l, clusters, params, HP);
                probs[prob_idx] = log_sum_exp(probs[prob_idx], msgs.a_msgs[a->id] + emission_ll);
            }
        }
        ++masked_l;
    }
    normalize_ll(probs, n_masked_ls, HP.K);
    return probs;
}
