#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include "max.hpp"
#include "seq_array.hpp"
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "math.hpp"
#include "util.hpp"


ViterbiBuffers::ViterbiBuffers(uint32_t n_cluster_ids, int L) :
    a_msgs(n_cluster_ids), new_a_msgs(L), new_b_msgs(L-1), path(2*L-1) {}

double get_cluster_emission_ll(
    Cluster* a, int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    if (xil == -1) {
        return 0.0;
    }

    if (a == nullptr) {
        if (clusters.emit_mode == EmitMode::soft) {
            return -std::log(static_cast<double>(HP.K));
        }

        int nkl = clusters.rs_by_emit[idx2d(l,xil,HP.K)].size();
        double ll = -delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K,
                                 clusters.rs[l].size(), params.mu_log_gamma[l]);
        if (clusters.emit_mode == EmitMode::noisy) {
            double c = (static_cast<double>(clusters.rs[l].size()) - HP.K*nkl) / (HP.K-1.0);
            double mu_y = params.mu_gamma[l] + nkl + c*params.mu_eps;
            double sigma2_y = params.sigma2_gamma[l] + c*c*params.sigma2_eps;
            ll += delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);
            return ll;
        }

        ll += delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nkl, params.mu_log_gamma[l]);
        return ll;
    }

    if (clusters.emit_mode == EmitMode::soft) {
        return delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0,
                           a->nk[xil], params.mu_log_gamma[l])
            - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n_obs, params.mu_log_gamma[l]);
    }
    if (clusters.emit_mode == EmitMode::noisy) {
        return xil == a->emission ? params.Eeps_log_match : params.Eeps_log_mismatch;
    }
    return xil == a->emission ? 0.0 : -std::numeric_limits<double>::infinity();
}

void get_viterbi_path(
    const SeqArray& x, int i, const std::unordered_map<int, int> *obs_ls,
    ViterbiBuffers& viterbi_bufs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<ViterbiMsg>& a_msgs = viterbi_bufs.a_msgs;
    a_msgs.resize(clusters.next_cluster_id);
    std::vector<ViterbiMsg>& new_a_msgs = viterbi_bufs.new_a_msgs;
    std::vector<ViterbiMsg>& new_b_msgs = viterbi_bufs.new_b_msgs;
    std::vector<Cluster*>& viterbi_path = viterbi_bufs.path;

    const std::vector<Cluster*> *matching_as = nullptr;
    for (int l = HP.L-1; l >= 0; --l) {
        int xil = get_xil(x, i, l, obs_ls);
        double new_a_ll = get_cluster_emission_ll(nullptr, xil, l, clusters, params, HP);

        matching_as = &clusters.get_matching_as(l, xil);
        if (l == HP.L-1) {
            new_a_msgs[l] = ViterbiMsg{new_a_ll, nullptr};
            for (Cluster *a : *matching_as) {
                a_msgs[a->id] = ViterbiMsg{get_cluster_emission_ll(a, xil, l, clusters, params, HP), nullptr};
            }
            continue;
        }

        // New b message.
        int nQl = clusters.qs[l].size();
        double mu_y = params.mu_alpha + nQl*params.mu_d[l];
        double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
        double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

        Cluster* best_a = nullptr;
        double best_a_ll = params.mu_log_alpha + new_a_msgs[l+1].ll;

        int xil1 = get_xil(x, i, l+1, obs_ls);
        for (Cluster *a : clusters.get_matching_as(l+1, xil1)) {
            double nCl = a->parents.size();
            double ll = params.mu_log_d[l] + std::log(nCl) + a_msgs[a->id].ll;
            if (ll > best_a_ll) {
                best_a = a;
                best_a_ll = ll;
            }
        }
        new_b_msgs[l] = ViterbiMsg{-elogy + best_a_ll, best_a};

        // a messages.
        new_a_msgs[l] = ViterbiMsg{new_a_ll + new_b_msgs[l].ll, nullptr};
        for (Cluster* a : *matching_as) {
            Cluster* best_b = nullptr;
            double nFl = a->children.size();
            double best_b_ll = std::log(nFl) + params.mu_log_d[l] + new_b_msgs[l].ll;
            for (Cluster* b : a->children) {
                Cluster* next_a = b->q_child;
                if (clusters.emit_mode == EmitMode::hard && xil1 != -1 && next_a->emission != xil1) {
                    continue;
                }
                double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n) + a_msgs[next_a->id].ll;
                if (ll > best_b_ll) {
                    best_b = b;
                    best_b_ll = ll;
                }
            }
            double emission_ll = get_cluster_emission_ll(a, xil, l, clusters, params, HP);
            a_msgs[a->id] = ViterbiMsg{emission_ll - std::log(static_cast<double>(a->n)) + best_b_ll, best_b};
        }
    }

    // Initial CRP probs.
    new_a_msgs[0].ll += params.mu_log_alpha;
    for (Cluster* a : *matching_as) {
        a_msgs[a->id].ll += std::log(static_cast<double>(a->n));
    }

    Cluster* a = nullptr;
    double a_ll = new_a_msgs[0].ll;
    for (Cluster* cand_a : *matching_as) {
        double ll = a_msgs[cand_a->id].ll;
        if (ll > a_ll) {
            a_ll = ll;
            a = cand_a;
        }
    }
    viterbi_path[0] = a;

    Cluster* b;
    for (int l = 0; l < HP.L-1; ++l) {
        b = a == nullptr ? new_a_msgs[l].next : a_msgs[a->id].next;
        viterbi_path[2 * l + 1] = b;

        a = b == nullptr ? new_b_msgs[l].next : b->q_child;
        viterbi_path[2 * l + 2] = a;
    }
}


int get_new_cluster_emission(
    int8_t xil, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    if (clusters.emit_mode == EmitMode::soft) {
        return -1;
    }
    if (xil == -1) {
        return clusters.cluster_mode(l);
    }
    if (clusters.emit_mode != EmitMode::noisy) {
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

void viterbi_add_seq(
    const SeqArray& x, int x_idx, int seq_idx,
    ViterbiBuffers& viterbi_bufs,
    Clusters& clusters, const Params& params, const HyperParams& HP
) {
    get_viterbi_path(x, x_idx, nullptr, viterbi_bufs, clusters, params, HP);
    std::vector<Cluster*>& viterbi_path = viterbi_bufs.path;

    Cluster* a = viterbi_path[0];
    Cluster* a_obj = a;
    if (a == nullptr) {
        int emission = get_new_cluster_emission(x(x_idx, 0), 0, clusters, params, HP);
        a_obj = clusters.create_empty_cluster(true, 0, emission);
    }
    clusters.cluster_add(a_obj, seq_idx, x(x_idx, 0));

    Cluster* b = nullptr;
    Cluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = viterbi_path[2 * l + 1];
        b_obj = (b == nullptr) ? clusters.create_empty_cluster(false, l, -1) : b;
        clusters.cluster_add(b_obj, seq_idx, -1);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = viterbi_path[2 * (l + 1)];
        a_obj = a;
        if (a == nullptr) {
            int emission = get_new_cluster_emission(x(x_idx, l+1), l+1, clusters, params, HP);
            a_obj = clusters.create_empty_cluster(true, l+1, emission);
        }
        clusters.cluster_add(a_obj, seq_idx, x(x_idx, l+1));
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
}

void max_step(const SeqArray& x, Clusters& clusters, const Params& params, const HyperParams& HP) {
    ViterbiBuffers viterbi_bufs(clusters.next_cluster_id, HP.L);
    for (int i = 0; i < HP.N; ++i) {
        for (int l = 0; l < HP.L; ++l) {
            clusters.cluster_remove(clusters.r_assign[idx2d(i, l, HP.L)], i, x(i, l));
            if (l == HP.L-1) {
                break;
            }
            clusters.cluster_remove(clusters.q_assign[idx2d(i, l, HP.L - 1)], i, -1);
        }
        viterbi_add_seq(x, i, i, viterbi_bufs, clusters, params, HP);
    }
}

void add_seqs(const SeqArray& x_new, Clusters& clusters, const Params& params, HyperParams& HP) {
    int old_N = HP.N;
    HP.N += x_new.N;
    clusters.r_assign.resize(HP.N * HP.L, nullptr);
    clusters.q_assign.resize(HP.N * (HP.L - 1), nullptr);

    ViterbiBuffers viterbi_bufs(clusters.next_cluster_id, HP.L);
    for (int i = 0; i < x_new.N; ++i) {
        viterbi_add_seq(x_new, i, old_N + i, viterbi_bufs, clusters, params, HP);
    }
}


void max_cluster_emissions(Clusters& clusters, const Params& params, const HyperParams& HP) {
    if (clusters.emit_mode != EmitMode::noisy) {
        throw std::runtime_error("Only need to maximize cluster emissions if noisy.");
    }
    for (int l = 0; l < HP.L; ++l) {
        for (Cluster* a : clusters.rs[l]) {
            int best_k = 0;
            double best_ll = -std::numeric_limits<double>::infinity();
            for (int k = 0; k < HP.K; ++k) {
                int nkl = clusters.rs_by_emit[idx2d(l,k,HP.K)].size() - (a->emission == k);
                double ll = delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l],
                                        1.0, nkl, params.mu_log_gamma[l])
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
