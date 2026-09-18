#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>
#include "fwd_bkwd.hpp"
#include "obs.hpp"
#include "seq_array.hpp"
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"
#include "math.hpp"
#include "model_array.hpp"
#include "util.hpp"


FwdBkwdMsgs::FwdBkwdMsgs(std::uint32_t n_cluster_ids, int L) : a(n_cluster_ids), new_a(L) {}

FwdBkwdBuffers::FwdBkwdBuffers(std::uint32_t n_cluster_ids, int L) :
    bkwd(n_cluster_ids, L), fwd(n_cluster_ids, L) {}


template <typename Obs>
void get_bkwd_msgs(
    const Obs& obs,
    FwdBkwdMsgs& msgs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<double>& a_msgs = msgs.a;
    std::vector<double>& new_a_msgs = msgs.new_a;
    for (int l = HP.L-1; l >= 0; --l) {
        double new_a_ll = obs.emission_ll(nullptr, l);

        const std::vector<Cluster*>& matching_as = obs.matching_as(l);
        if (l == HP.L-1) {
            new_a_msgs[l] = new_a_ll;
            for (Cluster* a : matching_as) {
                a_msgs[a->id] = obs.emission_ll(a, l);
            }
            continue;
        }

        // New b message.
        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        double new_b_ll = -elogy + params.mu_log_alpha + new_a_msgs[l+1];
        for (Cluster* a : obs.matching_as(l + 1)) {
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
                Cluster* next_a = b->q_child;
                if (!obs.matches(next_a, l + 1)) {
                    continue;
                }
                double ll = -log_na + delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n)
                    + a_msgs[next_a->id];
                a_ll = log_sum_exp(a_ll, ll);
            }
            a_msgs[a->id] = obs.emission_ll(a, l) + a_ll;
        }
    }
}

template <typename Obs>
void get_fwd_msgs(
    const Obs& obs,
    FwdBkwdMsgs& msgs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<double>& a_msgs = msgs.a;
    std::vector<double>& new_a_msgs = msgs.new_a;
    for (int l = 0; l < HP.L; ++l) {
        double new_a_ll = obs.emission_ll(nullptr, l);

        const std::vector<Cluster*>& matching_as = obs.matching_as(l);
        if (l == 0) {
            new_a_msgs[l] = new_a_ll + params.mu_log_alpha;
            for (Cluster* a : matching_as) {
                a_msgs[a->id] = obs.emission_ll(a, l)
                    + std::log(static_cast<double>(a->n));
            }
            continue;
        }

        // New b message.
        double new_b_ll = new_a_msgs[l-1];
        for (Cluster* a : obs.matching_as(l - 1)) {
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
                Cluster* prev_a = b->q_parent;
                if (!obs.matches(prev_a, l - 1)) {
                    continue;
                }
                double log_na = std::log(static_cast<double>(prev_a->n));
                double ab_ll = -log_na + delta_Elogx(params.mu_d[l-1], params.sigma2_d[l-1], -1.0, b->n);
                a_ll = log_sum_exp(a_ll, a_msgs[prev_a->id] + ab_ll);
            }
            a_msgs[a->id] = obs.emission_ll(a, l) + a_ll;
        }
    }
}

template <typename Obs>
void get_fwd_bkwd_msgs(
    const Obs& obs, FwdBkwdBuffers& fwd_bkwd_bufs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    FwdBkwdMsgs& msgs = fwd_bkwd_bufs.bkwd;
    FwdBkwdMsgs& fwd_msgs = fwd_bkwd_bufs.fwd;
    get_bkwd_msgs(obs, msgs, clusters, params, HP);
    get_fwd_msgs(obs, fwd_msgs, clusters, params, HP);

    for (int l = 0; l < HP.L; ++l) {
        msgs.new_a[l] += fwd_msgs.new_a[l] - obs.emission_ll(nullptr, l);
        for (Cluster* a : obs.matching_as(l)) {
            msgs.a[a->id] += fwd_msgs.a[a->id] - obs.emission_ll(a, l);
        }
    }
}

void fwd_bkwd(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    FwdBkwdBuffers& fwd_bkwd_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    int K = HP.n_emissions(0);
    FwdBkwdMsgs& msgs = fwd_bkwd_bufs.bkwd;
    ScalarObs obs{x, i, obs_ls, clusters, params, HP};
    get_fwd_bkwd_msgs(obs, fwd_bkwd_bufs, clusters, params, HP);

    int n_masked_ls = HP.L - obs_ls.size();
    std::fill(seq_probs.begin(), seq_probs.end(), -std::numeric_limits<double>::infinity());
    int masked_l = 0;
    for (int l = 0; l < HP.L; ++l) {
        if (obs_ls.contains(l)) { continue; }
        for (int k = 0; k < K; ++k) {
            std::size_t prob_idx = idx2d(masked_l, k, K);
            double emission_ll = get_new_cluster_emission_ll(k, l, clusters, params, HP);
            seq_probs[prob_idx] = log_sum_exp(seq_probs[prob_idx], msgs.new_a[l] + emission_ll);
        }
        for (Cluster* a : clusters.rs[l]) {
            std::size_t prob_idx = idx2d(masked_l, a->emission, K);
            seq_probs[prob_idx] = log_sum_exp(seq_probs[prob_idx], msgs.a[a->id]);
        }
        ++masked_l;
    }
    normalize_ll(seq_probs, n_masked_ls, K);
}

void fwd_bkwd_blocks(
    const ModelArray& ref, const BlockObs& obs,
    const std::unordered_map<int, int>& obs_ls,
    FwdBkwdBuffers& fwd_bkwd_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    get_fwd_bkwd_msgs(obs, fwd_bkwd_bufs, clusters, params, HP);
    FwdBkwdMsgs& msgs = fwd_bkwd_bufs.bkwd;

    std::fill(seq_probs.begin(), seq_probs.end(), -std::numeric_limits<double>::infinity());
    int masked_l = 0;
    std::vector<double> emission_state_ll;
    for (int l = 0; l < HP.L; ++l) {
        emission_state_ll.assign(
            HP.n_emissions(l), -std::numeric_limits<double>::infinity()
        );
        for (int emission = 0; emission < HP.n_emissions(l); ++emission) {
            if (!obs.is_emission_compatible(l, emission)) {
                continue;
            }
            emission_state_ll[emission] = msgs.new_a[l]
                + obs.new_emission_ll_at(l, emission)
                - obs.compatible_emission_ll[l];
        }
        for (Cluster* a : obs.matching_as(l)) {
            emission_state_ll[a->emission] = log_sum_exp(
                emission_state_ll[a->emission], msgs.a[a->id]
            );
        }

        for (int snp = ref.snp_start(l); snp < ref.snp_end(l); ++snp) {
            if (obs_ls.contains(snp)) {
                continue;
            }
            for (int emission = 0; emission < HP.n_emissions(l); ++emission) {
                int representative = ref.rep(l, emission);
                int allele = ref.snps(representative, snp);
                std::size_t prob_idx = idx2d(masked_l, allele, 2);
                seq_probs[prob_idx] = log_sum_exp(seq_probs[prob_idx], emission_state_ll[emission]);
            }
            ++masked_l;
        }
    }
    normalize_ll(seq_probs, masked_l, 2);
}
