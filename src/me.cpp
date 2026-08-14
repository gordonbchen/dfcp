#include <cstdint>
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <string_view>
#include <random>
#include <cstdlib>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "max.hpp"
#include "fwd_bkwd.hpp"
#include "expect.hpp"
#include "elbo.hpp"
#include "tree.hpp"
#include "json.hpp"
#include "util.hpp"


void train_dfcp(
    Clusters& clusters, Params& params, HyperParams& HP,

    bool block_init,
    bool pbwt_init, int pbwt_match_len, bool pbwt_match_curr,
    bool init_only,

    const std::vector<int8_t>& x_train,

    Json& json
) {
    // Init clusters.
    auto t0 = std::chrono::steady_clock::now();
    if (block_init) {
        clusters.block_init(x_train);
    }
    else if (pbwt_init) {
        if (HP.K != 2) { throw std::invalid_argument("pbwt_init only supported for K=2."); }
        clusters.pbwt_init(x_train, pbwt_match_len, pbwt_match_curr);
    }
    else {
        int N = HP.N;
        HP.N = 0;
        add_seqs(clusters, x_train.begin(), N, HP, params);
    }
    auto t1 = std::chrono::steady_clock::now();
    auto t_init = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cerr << "t_init=" << t_init << '\n';
    json.add("t_init", t_init);

    if (init_only) { return; }

    // Train.
    EarlyStopping early_stop{2, false, 1e-3};
    double elbo = 0.0;

    std::vector<Json> train_log;
    while (!early_stop.converged()) {
        auto t0 = std::chrono::steady_clock::now();
        max_step(clusters, x_train, HP, params);
        if (clusters.noisy) {
            max_cluster_emissions(clusters, HP, params);
        }
        auto t1 = std::chrono::steady_clock::now();
        expect_step(HP, params, clusters);
        auto t2 = std::chrono::steady_clock::now();
        elbo = calc_elbo(HP, params, clusters);
        auto t3 = std::chrono::steady_clock::now();
        early_stop.update(elbo);

        auto t_max = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        auto t_expect = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        auto t_elbo = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
        auto t_step = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t0).count();
        std::cerr << early_stop.step << ": elbo=" << elbo
            << " t_max=" << t_max << "ms t_expect=" << t_expect << "ms t_elbo=" << t_elbo
            << "ms t_step=" << t_step << "ms\n";
        train_log.emplace_back();
        train_log[train_log.size()-1].add("elbo", elbo)
            .add("t_max", t_max).add("t_expect", t_expect).add("t_elbo", t_elbo).add("t_step", t_step);
    }
    json.add("train_log", train_log);

    Json param_log;
    param_log.add("mu_alpha", params.mu_alpha).add("mu_gamma", params.mu_gamma).add("mu_d", params.mu_d)
        .add("alpha_eps", params.alpha_eps).add("beta_eps", params.beta_eps);
    json.add("params", param_log);
}


void eval_clusters(const Clusters& clusters, const HyperParams& HP, const std::vector<int8_t>& x_train, Json& json) {
    // Cluster stability IOU.
    auto t0 = std::chrono::steady_clock::now();
    double mean_iou = 0.0;
    double mean_emission_iou = 0.0;
    for (int l = 0; l < HP.L-1; ++l) {
        int n_intersect = 0;
        int n_union = 0;

        int n_emission_intersect = 0;
        int n_emission_union = 0;

        for (int i = 0; i < HP.N; ++i) {
            for (int j = i+1; j < HP.N; ++j) {
                int l_same = clusters.r_assign[idx2d(i,l,HP.L)] == clusters.r_assign[idx2d(j,l,HP.L)];
                int l1_same = clusters.r_assign[idx2d(i,l+1,HP.L)] == clusters.r_assign[idx2d(j,l+1,HP.L)];
                n_intersect += l_same && l1_same;
                n_union += l_same || l1_same;

                int l_emission_same = x_train[idx2d(i,l,HP.L)] == x_train[idx2d(j,l,HP.L)];
                int l1_emission_same = x_train[idx2d(i,l+1,HP.L)] == x_train[idx2d(j,l+1,HP.L)];
                n_emission_intersect += l_emission_same && l1_emission_same;
                n_emission_union += l_emission_same || l1_emission_same;
            }
        }
        mean_iou += (n_union == 0) ? 1.0 : static_cast<double>(n_intersect) / n_union;
        mean_emission_iou += (n_emission_union == 0) ? 1.0
            : static_cast<double>(n_emission_intersect) / n_emission_union;
    }
    mean_iou /= HP.L-1;
    mean_emission_iou /= HP.L-1;
    auto t1 = std::chrono::steady_clock::now();
    auto t_iou = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cerr << "mean_iou=" << mean_iou << " mean_emission_iou=" << mean_emission_iou
        << " t_iou=" << t_iou << "ms\n";
    json.add("mean_iou", mean_iou).add("mean_emission_iou", mean_emission_iou).add("t_iou", t_iou);

    // Average # of clusters.
    double mean_clusters = static_cast<double>(clusters.nR) / HP.L;
    std::cerr << "mean_clusters=" << mean_clusters << '\n';
    json.add("mean_clusters", mean_clusters);

    // Purity.
    double cluster_purity = 1.0;
    if (clusters.soft) {
        cluster_purity = 0.0;
        for (int l = 0; l < HP.L; ++l) {
            for (Cluster *c : clusters.rs[l]) {
                cluster_purity += c->mode().count;
            }
        }
        cluster_purity /= HP.N*HP.L;
    }
    else if (clusters.noisy) {
        cluster_purity = clusters.n_matches / static_cast<double>(clusters.n_obs);
    }
    std::cerr << "cluster_purity=" << cluster_purity << '\n';
    json.add("cluster_purity", cluster_purity);
}


void tree_eval(
    const Clusters& clusters, const HyperParams& HP,

    const std::vector<int8_t>& x_train, const std::vector<size_t>& train_idxs, int n_val_seqs,

    char *tree_fname, char *variant_pos_fname, int variant_start_pos,
    double clade_beta,
    char *tree_vis_fname,

    Json& json
) {
    if ((tree_fname == nullptr) || (variant_pos_fname == nullptr) || (variant_start_pos < 0)) {
        throw std::invalid_argument("Tree eval requires variant position file and variant start pos.");
    }
    std::vector<int> variant_pos{parse_pos_file_idx(variant_pos_fname, variant_start_pos, HP.L)};
    auto [coal_trees, recomb_pos] = parse_tree_file(tree_fname);
    std::vector<int> tree_idxs{get_tree_idxs(variant_pos, recomb_pos)};

    std::unordered_set<size_t> train_idxs_set(train_idxs.begin(), train_idxs.end());

    int excess_parsimony = 0;
    int emission_excess_parsimony = 0;
    std::chrono::steady_clock::duration t_parsimony_dur{};

    double weighted_clade_iou_sum = 0.0;
    double clade_weight_sum = 0.0;
    std::vector<std::vector<double>> dfcp_clade_heights(HP.L);
    double emission_weighted_clade_iou_sum = 0.0;
    double emission_clade_weight_sum = 0.0;
    std::vector<std::vector<double>> emission_clade_heights(HP.L);
    std::vector<double> coal_root_heights(HP.L);
    std::chrono::steady_clock::duration t_clade_iou_dur{};

    for (int l = 0; l < HP.L; ++l) {
        const std::unordered_map<int, CoalNode>& coal_tree{coal_trees[tree_idxs[l]]};
        coal_root_heights[l] = calc_node_height(coal_tree, -1);
        dfcp_clade_heights[l].reserve(clusters.rs[l].size());
        emission_clade_heights[l].reserve(HP.K);

        // Parismony.
        auto t_parsimony0 = std::chrono::steady_clock::now();

        std::unordered_map<Cluster*, int> cluster_idxs;
        cluster_idxs.reserve(clusters.rs[l].size());
        for (Cluster* c : clusters.rs[l]) {
            cluster_idxs.emplace(c, cluster_idxs.size());
        }
        std::vector<int> cluster_assign(HP.N + n_val_seqs, -1);
        for (int i = 0; i < HP.N; ++i) {
            cluster_assign[train_idxs[i]] = cluster_idxs.at(clusters.r_assign[idx2d(i, l, HP.L)]);
        }
        excess_parsimony += calc_excess_parsimony(coal_tree, cluster_assign, clusters.rs[l].size(),
                                                  clusters.rs[l].size(), train_idxs_set);

        std::vector<int> emission_counts(HP.K, 0);
        std::vector<int> emission_clusters(HP.N + n_val_seqs, -1);
        for (int i = 0; i < HP.N; ++i) {
            emission_clusters[train_idxs[i]] = x_train[idx2d(i,l,HP.L)];
            ++emission_counts[x_train[idx2d(i,l,HP.L)]];
        }
        int n_obs_emissions = std::count_if(emission_counts.begin(), emission_counts.end(),
                                            [](int count) { return count > 0; });
        emission_excess_parsimony += calc_excess_parsimony(coal_tree, emission_clusters, n_obs_emissions,
                                                           HP.K, train_idxs_set);

        t_parsimony_dur += std::chrono::steady_clock::now() - t_parsimony0;

        // Importance-weighted clade iou.
        auto t_clade_iou0 = std::chrono::steady_clock::now();
        for (const auto& [c, cluster_idx] : cluster_idxs) {
            auto [max_clade_iou, clade_root] = calc_max_clade_iou(coal_tree, cluster_assign,
                                                                  cluster_idx, c->n, train_idxs_set);
            dfcp_clade_heights[l].emplace_back(calc_node_height(coal_tree, clade_root));

            double z = static_cast<double>(c->n - 1) / (HP.N - 1);
            double weight = std::pow(z * (1.0 - z), clade_beta - 1.0);
            weighted_clade_iou_sum += weight * max_clade_iou;
            clade_weight_sum += weight;
        }

        for (int k = 0; k < HP.K; ++k) {
            if (emission_counts[k] == 0) {
                continue;
            }
            auto [max_clade_iou, clade_root] = calc_max_clade_iou(coal_tree, emission_clusters,
                                                                  k, emission_counts[k], train_idxs_set);
            emission_clade_heights[l].emplace_back(calc_node_height(coal_tree, clade_root));

            double z = static_cast<double>(emission_counts[k] - 1) / (HP.N - 1);
            double weight = std::pow(z * (1.0 - z), clade_beta - 1.0);
            emission_weighted_clade_iou_sum += weight * max_clade_iou;
            emission_clade_weight_sum += weight;
        }
        t_clade_iou_dur += std::chrono::steady_clock::now() - t_clade_iou0;

        // Tree viz.
        if ((tree_vis_fname != nullptr) && (l < 16)) {
            tree_to_dot(tree_vis_fname, coal_tree, cluster_assign, emission_clusters, train_idxs_set, l, 16);
        }
    }
    double mean_excess_parsimony = static_cast<double>(excess_parsimony) / HP.L;
    double mean_emission_excess_parsimony = static_cast<double>(emission_excess_parsimony) / HP.L;

    double clade_iou = clade_weight_sum == 0.0 ? -1.0 : weighted_clade_iou_sum / clade_weight_sum;
    if (clade_weight_sum == 0.0) { std::cerr << "clade_iou undefined, no nontrivial clusters.\n"; }
    double emission_clade_iou = emission_clade_weight_sum == 0.0 ?
        -1.0 : emission_weighted_clade_iou_sum / emission_clade_weight_sum;
    if (emission_clade_weight_sum == 0.0) {std::cerr << "emission_clade_iou undefined, no nontrivial clusters.\n";}

    auto t_parsimony = std::chrono::duration_cast<std::chrono::milliseconds>(t_parsimony_dur).count();
    auto t_clade_iou = std::chrono::duration_cast<std::chrono::milliseconds>(t_clade_iou_dur).count();
    std::cerr << "mean_excess_parsimony=" << mean_excess_parsimony
        << " mean_emission_excess_parsimony=" << mean_emission_excess_parsimony
        << " t_parsimony=" << t_parsimony << "ms\n"
        << "clade_iou=" << clade_iou << " emission_clade_iou=" << emission_clade_iou
        << " clade_beta=" << clade_beta << " t_clade_iou=" << t_clade_iou << "ms\n";
    json.add("mean_excess_parsimony", mean_excess_parsimony)
        .add("mean_emission_excess_parsimony", mean_emission_excess_parsimony)
        .add("t_parsimony", t_parsimony)
        .add("clade_iou", clade_iou).add("emission_clade_iou", emission_clade_iou)
        .add("clade_beta", clade_beta).add("t_clade_iou", t_clade_iou)
        .add("dfcp_clade_heights", dfcp_clade_heights)
        .add("emission_clade_heights", emission_clade_heights)
        .add("coal_root_heights", coal_root_heights);
}


int get_max_prob_emission(std::vector<double>::const_iterator probs, int K) {
    int best_k = 0;
    double best_prob = probs[0];
    for (int k = 1; k < K; ++k) {
        if (probs[k] > best_prob) {
            best_k = k;
            best_prob = probs[k];
        }
    }
    return best_k;
}

struct R2Calc {
    int n = 0;
    double mean_x = 0.0;
    double mean_y = 0.0;
    double var_x = 0.0;
    double var_y = 0.0;
    double cov = 0.0;

    void update(double x, double y) {
        ++n;
        double delta_x = x - mean_x;
        double delta_y = y - mean_y;
        mean_x += delta_x / n;
        mean_y += delta_y / n;
        var_x += delta_x * (x - mean_x);
        var_y += delta_y * (y - mean_y);
        cov += delta_x * (y - mean_y);
    }

    double r2() const {
        if ((var_x <= 0.0) || (var_y <= 0.0)) {
            return -1.0;
        }
        return std::clamp((cov * cov) / (var_x * var_y), 0.0, 1.0);
    }
};

void impute(
    const Clusters& clusters, const Params& params, const HyperParams& HP,

    const std::vector<int8_t>& x_train,
    const std::vector<int8_t>& x_val, int n_val_seqs,
    const std::vector<int8_t>& x_val_true, const std::vector<int>& masked_ls, int n_masked_ls,

    Json& json
) {
    std::vector<size_t> emission_counts{count_emissions(x_train, HP.N, HP.L, HP.K)};
    std::vector<int8_t> modes{get_emission_modes(emission_counts, HP.L, HP.K)};

    int n_viterbi_correct = 0;
    int n_fwd_bkwd_correct = 0;
    int n_mode_correct = 0;

    int n_r2 = (HP.K == 2) ? n_masked_ls : 0;
    std::vector<R2Calc> fwd_bkwd_correlations(n_r2);
    std::vector<R2Calc> viterbi_correlations(n_r2);
    auto [minor_alleles, minor_allele_counts] = count_minor_alleles(emission_counts, masked_ls, n_masked_ls, HP.K);

    std::chrono::steady_clock::duration viterbi_impute_dur{};
    std::chrono::steady_clock::duration fwd_bkwd_impute_dur{};
    for (int i = 0; i < n_val_seqs; ++i) {
        auto t0 = std::chrono::steady_clock::now();
        std::vector<double> fwd_bkwd_probs{fwd_bkwd(x_val.begin() + i*HP.L, masked_ls, clusters, params, HP)};
        fwd_bkwd_impute_dur += std::chrono::steady_clock::now() - t0;

        t0 = std::chrono::steady_clock::now();
        std::vector<Cluster*> viterbi_clusters{get_viterbi_clusters(clusters, x_val.begin() + i*HP.L, HP, params)};
        std::vector<double> viterbi_probs;
        viterbi_probs.reserve(n_masked_ls * HP.K);
        for (int il = 0; il < n_masked_ls; ++il) {
            int l = masked_ls[il];
            for (int k = 0; k < HP.K; ++k) {
                double p = get_cluster_emission_ll(viterbi_clusters[2*l], k, l, clusters, params, HP);
                viterbi_probs.emplace_back(p);
            }
        }
        normalize_ll(viterbi_probs, n_masked_ls, HP.K);
        viterbi_impute_dur += std::chrono::steady_clock::now() - t0;

        for (int il = 0; il < n_masked_ls; ++il) {
            int l = masked_ls[il];
            int x_true = x_val_true[idx2d(i,il,n_masked_ls)];

            // Accuracy.
            if (x_true == get_max_prob_emission(viterbi_probs.begin() + il*HP.K, HP.K)) {
                ++n_viterbi_correct;
            }
            if (x_true == get_max_prob_emission(fwd_bkwd_probs.begin() + il*HP.K, HP.K)) {
                ++n_fwd_bkwd_correct;
            }
            if (x_true == modes[l]) {
                ++n_mode_correct;
            }

            // r^2 b/t minor allele and p(minor allele).
            if (HP.K == 2) {
                int minor_allele = minor_alleles[il];
                int is_minor_allele = x_true == minor_allele;
                fwd_bkwd_correlations[il].update(is_minor_allele, fwd_bkwd_probs[idx2d(il,minor_allele, HP.K)]);
                viterbi_correlations[il].update(is_minor_allele, viterbi_probs[idx2d(il,minor_allele, HP.K)]);
            }
        }
    }
    double n_masked_alleles = n_val_seqs * n_masked_ls;
    double viterbi_impute_acc = n_viterbi_correct / n_masked_alleles;
    double mode_impute_acc = n_mode_correct / n_masked_alleles;
    double fwd_bkwd_impute_acc = n_fwd_bkwd_correct / n_masked_alleles;

    std::vector<double> fwd_bkwd_r2s(n_r2, -1.0);
    std::vector<double> viterbi_r2s(n_r2, -1.0);
    double fwd_bkwd_r2_mean = 0.0;
    double viterbi_r2_mean = 0.0;
    int n_fwd_bkwd_r2_locs = 0;
    int n_viterbi_r2_locs = 0;
    if (HP.K == 2) {
        for (int il = 0; il < n_masked_ls; ++il) {
            fwd_bkwd_r2s[il] = fwd_bkwd_correlations[il].r2();
            viterbi_r2s[il] = viterbi_correlations[il].r2();
            if (fwd_bkwd_r2s[il] >= 0.0) {
                fwd_bkwd_r2_mean += fwd_bkwd_r2s[il];
                ++n_fwd_bkwd_r2_locs;
            }
            if (viterbi_r2s[il] >= 0.0) {
                viterbi_r2_mean += viterbi_r2s[il];
                ++n_viterbi_r2_locs;
            }
        }
        fwd_bkwd_r2_mean = (n_fwd_bkwd_r2_locs == 0) ? -1.0 : fwd_bkwd_r2_mean / n_fwd_bkwd_r2_locs;
        viterbi_r2_mean = (n_viterbi_r2_locs == 0) ? -1.0 : viterbi_r2_mean / n_viterbi_r2_locs;
        std::cerr << "fwd_bkwd_r2_mean=" << fwd_bkwd_r2_mean << " viterbi_r2_mean=" << viterbi_r2_mean << '\n';
        json.add("fwd_bkwd_r2s", fwd_bkwd_r2s).add("viterbi_r2s", viterbi_r2s)
            .add("fwd_bkwd_r2_mean", fwd_bkwd_r2_mean).add("viterbi_r2_mean", viterbi_r2_mean)
            .add("minor_allele_counts", minor_allele_counts);
    }

    auto t_viterbi_impute = std::chrono::duration_cast<std::chrono::milliseconds>(viterbi_impute_dur).count();
    auto t_fwd_bkwd_impute = std::chrono::duration_cast<std::chrono::milliseconds>(fwd_bkwd_impute_dur).count();
    std::cerr << "viterbi_impute_acc=" << viterbi_impute_acc << " t_viterbi_impute=" << t_viterbi_impute << "ms\n"
        << "fwd_bkwd_impute_acc=" << fwd_bkwd_impute_acc << " t_fwd_bkwd_impute=" << t_fwd_bkwd_impute << "ms\n"
        << "mode_impute_acc=" << mode_impute_acc << '\n';
    json.add("viterbi_impute_acc", viterbi_impute_acc).add("t_viterbi_impute", t_viterbi_impute)
        .add("fwd_bkwd_impute_acc", fwd_bkwd_impute_acc).add("t_fwd_bkwd_impute", t_fwd_bkwd_impute)
        .add("mode_impute_acc", mode_impute_acc);
}


int main(int argc, char *argv[]) {
    Json json;

    // Read seq file.
    if (argc < 2) { throw std::invalid_argument("Requires sequence file."); }
    std::cerr << "seq_file=" << argv[1] << '\n';
    json.add("seq_file", argv[1]);
    std::ifstream seq_file(argv[1]);
    if (!seq_file.is_open()) { throw std::runtime_error("Failed to open sequence file."); };

    int N = 0;
    int L = 0;
    std::vector<int8_t> x_raw;
    std::string line;
    while (std::getline(seq_file, line)) {
        if (N == 0) {
            L = line.length();
        }
        for (char c : line) {
            if (c < '0' || c > '9') { throw std::runtime_error("Invalid allele char."); };
            x_raw.push_back(c - '0');
        }
        ++N;
    }
    int K = *std::max_element(x_raw.begin(), x_raw.end()) + 1;
    HyperParams HP{.N=N, .L=L, .K=K};

    // Parse optional args.
    int seed = 0;
    double val = 0.0;
    double mask = 0.0;

    char *tree_fname = nullptr;
    char *variant_pos_fname = nullptr;
    int variant_start_pos = -1;
    char *tree_vis_fname = nullptr;
    double clade_beta = 2.0;

    bool noisy = false;
    bool soft = false;
    bool block_init = false;
    bool pbwt_init = false;
    int pbwt_match_len = 5;
    bool pbwt_match_curr = true;
    bool init_only = false;

    int i = 2;
    while (i < argc) {
        if (i+1 >= argc) { throw std::invalid_argument("Arg has no value."); };

        std::string_view arg{argv[i]};
        if (arg == "--tau_1") { HP.tau_1 = parse_double(argv[i+1]); }
        else if (arg == "--tau_2") { HP.tau_2 = parse_double(argv[i+1]); }
        else if (arg == "--v_1") { HP.v_1 = parse_double(argv[i+1]); }
        else if (arg == "--v_2") { HP.v_2 = parse_double(argv[i+1]); }
        else if (arg == "--phi_1") { HP.phi_1 = parse_double(argv[i+1]); }
        else if (arg == "--phi_2") { HP.phi_2 = parse_double(argv[i+1]); }

        else if (arg == "--noisy") { noisy = parse_int(argv[i+1]) == 1; }
        else if (arg == "--lambda_1") { HP.lambda_1 = parse_double(argv[i+1]); }
        else if (arg == "--lambda_2") { HP.lambda_2 = parse_double(argv[i+1]); }

        else if (arg == "--seed") { seed = parse_int(argv[i+1]); }
        else if (arg == "--val") {
            val = parse_double(argv[i+1]);
            if ((val < 0.0) || (val >= 1.0)) { throw std::invalid_argument("val must be in [0, 1)."); }
        }
        else if (arg == "--mask") {
            mask = parse_double(argv[i+1]);
            if ((mask < 0.0) || (mask > 1.0)) { throw std::invalid_argument("mask must be in [0, 1]."); }
        }

        else if (arg == "--tree") { tree_fname = argv[i+1]; }
        else if (arg == "--variant_pos_fname") { variant_pos_fname = argv[i+1]; }
        else if (arg == "--variant_start_pos") { variant_start_pos = parse_int(argv[i+1]); }
        else if (arg == "--tree_vis") { tree_vis_fname = argv[i+1]; }
        else if (arg == "--clade_beta") {
            clade_beta = parse_double(argv[i+1]);
            if (clade_beta < 1.0) { throw std::invalid_argument("clade_beta must be at least 1."); }
        }

        else if (arg == "--soft") { soft = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--block_init") { block_init = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--pbwt_init") { pbwt_init = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--pbwt_match_len") { pbwt_match_len = parse_int(argv[i+1]); }
        else if (arg == "--pbwt_match_curr") { pbwt_match_curr = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--init_only") { init_only = (parse_int(argv[i+1]) == 1); }

        else { throw std::invalid_argument("Arg not recognized."); }
        i += 2;
    }
    if (noisy && soft) { throw std::invalid_argument("noisy is only for hard dfcp."); }
    if (block_init && pbwt_init) { throw std::invalid_argument("cannot do block and pbwt init."); }
    std::cerr << HP << '\n';

    // Split val for imputation.
    if ((val > 0.0) != (mask > 0.0)) {
        throw std::invalid_argument("If imputation val frac > 0, need mask frac > 0.");
    };
    std::vector<int8_t> x_train;
    std::vector<size_t> train_idxs;

    int n_val_seqs = 0;
    std::vector<int8_t> x_val;

    int n_masked_ls = 0;
    std::vector<int> masked_ls;
    std::vector<int8_t> x_val_true;

    if (val > 0.0) {
        std::random_device rd;
        seed = (seed == 0) ? rd() : seed;
        std::cerr << "seed=" << seed << '\n';
        json.add("seed", seed);
        std::mt19937 gen(seed);

        n_val_seqs = static_cast<int>(std::lround(N * val));
        if (n_val_seqs == 0) { throw std::invalid_argument("no val seqs."); }
        if (n_val_seqs == N) { throw std::invalid_argument("all val seqs."); }
        std::vector<int> seq_idxs(N);
        std::iota(seq_idxs.begin(), seq_idxs.end(), 0);
        std::vector<int> val_idxs{sample_without_replacement(std::move(seq_idxs), n_val_seqs, gen)};
        std::vector<bool> is_val(N, false);
        for (int i : val_idxs) {
            is_val[i] = true;
        }

        int n_train_seqs = N - n_val_seqs;
        x_train.reserve(n_train_seqs * HP.L);
        train_idxs.reserve(n_train_seqs);
        x_val.reserve(n_val_seqs * HP.L);

        for (int i = 0; i < N; ++i) {
            auto line = x_raw.begin() + i*HP.L;
            if (is_val[i]) {
                x_val.insert(x_val.end(), line, line+HP.L);
            }
            else {
                x_train.insert(x_train.end(), line, line+HP.L);
                train_idxs.emplace_back(i);
            }
        }

        HP.N = n_train_seqs;
        std::vector<size_t> train_emission_counts{count_emissions(x_train, HP.N, HP.L, HP.K)};
        std::vector<int> mask_eligible_ls;
        mask_eligible_ls.reserve(HP.L);
        for (int l = 0; l < HP.L; ++l) {
            auto l_counts = train_emission_counts.begin() + l*HP.K;
            int n_observed_alleles = std::count_if(
                l_counts, l_counts + HP.K, [](size_t count) { return count > 0; }
            );
            if (n_observed_alleles >= 2) {
                mask_eligible_ls.emplace_back(l);
            }
        }
        n_masked_ls = static_cast<int>(std::lround(HP.L * mask));
        if (n_masked_ls == 0) { throw std::invalid_argument("no masked ls."); }
        if (n_masked_ls > static_cast<int>(mask_eligible_ls.size())) {
            throw std::invalid_argument("not enough training-polymorphic loci for requested mask fraction.");
        }
        masked_ls = sample_without_replacement(std::move(mask_eligible_ls), n_masked_ls, gen);

        x_val_true.reserve(n_val_seqs * n_masked_ls);
        for (int i = 0; i < n_val_seqs; ++i) {
            for (int l : masked_ls) {
                x_val_true.emplace_back(x_val[idx2d(i, l, HP.L)]);
                x_val[idx2d(i, l, HP.L)] = -1;
            }
        }
        std::cerr << HP << "\nn_val_seqs=" << n_val_seqs << " n_masked_ls=" << n_masked_ls << '\n';
        json.add("n_val_seqs", n_val_seqs).add("n_masked_ls", n_masked_ls);
    }
    else {
        x_train = std::move(x_raw);
        for (int i = 0; i < HP.N; ++i) {
            train_idxs.emplace_back(i);
        }
    }

    // Init params and clusters.
    Params params{HP};
    Clusters clusters{HP, soft, noisy};

    train_dfcp(
        clusters, params, HP,

        block_init,
        pbwt_init, pbwt_match_len, pbwt_match_curr,
        init_only,

        x_train,

        json
    );

    eval_clusters(clusters, HP, x_train, json);
    if (tree_fname != nullptr) {
        tree_eval(
            clusters, HP,

            x_train, train_idxs, n_val_seqs,

            tree_fname, variant_pos_fname, variant_start_pos,
            clade_beta,
            tree_vis_fname,

            json
        );
    }

    if (n_val_seqs > 0) {
        impute(
            clusters, params, HP,

            x_train,
            x_val, n_val_seqs,
            x_val_true, masked_ls, n_masked_ls,

            json
        );
    }
    std::cout << json.str() << '\n';
}
