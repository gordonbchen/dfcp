#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <omp.h>
#include <unordered_map>
#include <vector>
#include "obs.hpp"
#include "max.hpp"
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "math.hpp"
#include "model_array.hpp"
#include "seq_array.hpp"
#include "util.hpp"


ViterbiBuffers::ViterbiBuffers(uint32_t n_cluster_ids, int L) :
    a_msgs(n_cluster_ids), new_a_msgs(L), new_b_msgs(L-1), path(2*L-1) {}


template <typename Obs>
void get_viterbi_path_impl(
    const Obs& obs,
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
        double new_a_ll = obs.emission_ll(nullptr, l);

        matching_as = &obs.matching_as(l);
        if (l == HP.L-1) {
            new_a_msgs[l] = ViterbiMsg{new_a_ll, nullptr};
            for (Cluster *a : *matching_as) {
                a_msgs[a->id] = ViterbiMsg{obs.emission_ll(a, l), nullptr};
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

        for (Cluster *a : obs.matching_as(l + 1)) {
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
                if (!obs.matches(next_a, l + 1)) {
                    continue;
                }
                double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->n) + a_msgs[next_a->id].ll;
                if (ll > best_b_ll) {
                    best_b = b;
                    best_b_ll = ll;
                }
            }
            double emission_ll = obs.emission_ll(a, l);
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

void get_viterbi_path(
    const ScalarObs& obs,
    ViterbiBuffers& viterbi_bufs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    get_viterbi_path_impl(obs, viterbi_bufs, clusters, params, HP);
}

void get_viterbi_path(
    const BlockObs& obs, ViterbiBuffers& viterbi_bufs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    get_viterbi_path_impl(obs, viterbi_bufs, clusters, params, HP);
}


void viterbi_add_path(
    const ModelArray& x, int x_idx, int seq_idx, Cluster* const* viterbi_path,
    Clusters& clusters, const HyperParams& HP
) {
    Cluster* a = viterbi_path[0];
    Cluster* a_obj = a;
    if (a == nullptr) {
        a_obj = clusters.create_empty_cluster(true, 0, x(x_idx, 0));
    }
    clusters.cluster_add(a_obj, seq_idx);

    Cluster* b = nullptr;
    Cluster* b_obj = nullptr;
    for (int l = 0; l < HP.L-1; ++l) {
        b = viterbi_path[2 * l + 1];
        b_obj = (b == nullptr) ? clusters.create_empty_cluster(false, l, -1) : b;
        clusters.cluster_add(b_obj, seq_idx);
        if (a == nullptr || b == nullptr) {
            a_obj->add_child(b_obj);
        }

        a = viterbi_path[2 * (l + 1)];
        a_obj = a;
        if (a == nullptr) {
            a_obj = clusters.create_empty_cluster(true, l+1, x(x_idx, l+1));
        }
        clusters.cluster_add(a_obj, seq_idx);
        if (a == nullptr || b == nullptr) {
            b_obj->add_child(a_obj);
        }
    }
}

void remove_seq(int i, Clusters& clusters, const HyperParams& HP) {
    for (int l = 0; l < HP.L; ++l) {
        clusters.cluster_remove(clusters.r_assign[idx2d(i, l, HP.L)], i);
        if (l < HP.L-1) {
            clusters.cluster_remove(clusters.q_assign[idx2d(i, l, HP.L - 1)], i);
        }
    }
}

void max_step(
    const ModelArray& x, Clusters& clusters, const Params& params, const HyperParams& HP,
    int batch_size
) {
    if (batch_size == 1) {
        ViterbiBuffers viterbi_bufs(clusters.next_cluster_id, HP.L);
        for (int i = 0; i < HP.N; ++i) {
            remove_seq(i, clusters, HP);
            ScalarObs obs{x, i, clusters, params, HP};
            get_viterbi_path(obs, viterbi_bufs, clusters, params, HP);
            viterbi_add_path(x, i, i, viterbi_bufs.path.data(), clusters, HP);
        }
        return;
    }

    int n_threads = std::min(omp_get_max_threads(), batch_size);
    std::vector<ViterbiBuffers> thread_bufs;
    thread_bufs.reserve(n_threads);
    for (int thread = 0; thread < n_threads; ++thread) {
        thread_bufs.emplace_back(clusters.next_cluster_id, HP.L);
    }

    std::size_t path_size = static_cast<std::size_t>(2*HP.L - 1);
    std::vector<Cluster*> paths(static_cast<std::size_t>(batch_size) * path_size);

    for (int i = 0; i < HP.N; i += batch_size) {
        int size = std::min(batch_size, HP.N - i);
        for (int j = 0; j < size; ++j) {
            remove_seq(i + j, clusters, HP);
        }

        #pragma omp parallel for num_threads(n_threads)
        for (int j = 0; j < size; ++j) {
            ViterbiBuffers& bufs = thread_bufs[omp_get_thread_num()];
            ScalarObs obs{x, i + j, clusters, params, HP};
            get_viterbi_path(obs, bufs, clusters, params, HP);
            std::copy(bufs.path.begin(), bufs.path.end(), paths.begin() + j*path_size);
        }

        for (int j = 0; j < size; ++j) {
            viterbi_add_path(x, i+j, i+j, paths.data() + j*path_size, clusters, HP);
        }
    }
}

void add_seqs(const ModelArray& x_new, Clusters& clusters, const Params& params, HyperParams& HP) {
    int old_N = HP.N;
    HP.N += x_new.snps.N;
    clusters.r_assign.resize(HP.N * HP.L, nullptr);
    clusters.q_assign.resize(HP.N * (HP.L - 1), nullptr);

    ViterbiBuffers viterbi_bufs(clusters.next_cluster_id, HP.L);
    for (int i = 0; i < x_new.snps.N; ++i) {
        ScalarObs obs{x_new, i, clusters, params, HP};
        get_viterbi_path(obs, viterbi_bufs, clusters, params, HP);
        viterbi_add_path(x_new, i, old_N + i, viterbi_bufs.path.data(), clusters, HP);
    }
}


void get_viterbi_impute_probs(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    ViterbiBuffers& viterbi_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    int K = HP.n_emissions(0);
    ScalarObs obs{x, i, obs_ls, clusters, params, HP};
    get_viterbi_path(obs, viterbi_bufs, clusters, params, HP);
    int n_masked_ls = HP.L - obs_ls.size();
    int masked_l = 0;
    for (int l = 0; l < HP.L; ++l) {
        if (obs_ls.contains(l)) { continue; }
        Cluster* a = viterbi_bufs.path[2 * l];
        if (a == nullptr) {
            for (int k = 0; k < K; ++k) {
                seq_probs[idx2d(masked_l, k, K)] = get_new_cluster_emission_ll(k, l, clusters, params, HP);
            }
        }
        else {
            std::fill_n(seq_probs.begin() + idx2d(masked_l, 0, K), K,
                        -std::numeric_limits<double>::infinity());
            seq_probs[idx2d(masked_l, a->emission, K)] = 0.0;
        }
        ++masked_l;
    }
    normalize_ll(seq_probs, n_masked_ls, K);
}

void get_blocked_viterbi_impute_probs(
    const ModelArray& ref, const BlockObs& obs, const std::unordered_map<int, int>& obs_ls,
    ViterbiBuffers& viterbi_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    get_viterbi_path(obs, viterbi_bufs, clusters, params, HP);
    int masked_snp = 0;
    for (int l = 0; l < HP.L; ++l) {
        Cluster* a = viterbi_bufs.path[2 * l];
        int emission = a == nullptr ? -1 : a->emission;
        if (emission == -1) {
            double best_ll = -std::numeric_limits<double>::infinity();
            for (int candidate = 0; candidate < HP.n_emissions(l); ++candidate) {
                if (!obs.is_emission_compatible(l, candidate)) {
                    continue;
                }
                double ll = obs.new_emission_ll_at(l, candidate);
                if (ll > best_ll) {
                    emission = candidate;
                    best_ll = ll;
                }
            }
        }
        int representative = ref.rep(l, emission);
        for (int snp = ref.snp_start(l); snp < ref.snp_end(l); ++snp) {
            if (obs_ls.contains(snp)) {
                continue;
            }
            int allele = ref.snps(representative, snp);
            seq_probs[idx2d(masked_snp, 0, 2)] = allele == 0 ? 1.0 : 0.0;
            seq_probs[idx2d(masked_snp, 1, 2)] = allele == 1 ? 1.0 : 0.0;
            ++masked_snp;
        }
    }
}
