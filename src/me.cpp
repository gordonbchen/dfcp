#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <unordered_map>
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

    std::vector<int8_t>& x_train,

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


void dfcp(
    HyperParams& HP,

    bool soft, bool noisy,

    bool block_init,
    bool pbwt_init, int pbwt_match_len, bool pbwt_match_curr,
    bool init_only,

    std::vector<int8_t>& x_train, int n_train_seqs, std::vector<size_t>& train_idxs,
    std::vector<int8_t>& x_val, int n_val_seqs, std::vector<size_t>& val_idxs,
    std::vector<SparseX>& x_val_true, int n_masked_alleles,

    char *tree_fname, char *variant_pos_fname, int variant_start_pos,
    char *tree_vis_fname, double clade_beta,

    Json& json
) {
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

    // Impute.
    if (n_val_seqs > 0) {
        auto t0 = std::chrono::steady_clock::now();
        add_seqs(clusters, x_val.begin(), n_val_seqs, HP, params);
        auto t1 = std::chrono::steady_clock::now();
        auto t_impute = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        std::vector<int8_t> modes{count_modes(x_train, n_train_seqs, HP.L, HP.K)};

        int n_dfcp_correct = 0;
        int n_mode_correct = 0;
        for (SparseX& s : x_val_true) {
            if (s.x == clusters.r_assign[idx2d(n_train_seqs+s.i, s.l, HP.L)]->get_imputed_emission(clusters.soft)) {
                ++n_dfcp_correct;
            }
            if (s.x == modes[s.l]) {
                ++n_mode_correct;
            }
        }
        double dfcp_impute_acc = static_cast<double>(n_dfcp_correct) / n_masked_alleles;
        double mode_impute_acc = static_cast<double>(n_mode_correct) / n_masked_alleles;

        std::cerr << "dfcp_impute_acc=" << dfcp_impute_acc << " t_impute=" << t_impute
            << "ms mode_impute_acc=" << mode_impute_acc << '\n';
        json.add("dfcp_impute_acc", dfcp_impute_acc).add("t_impute", t_impute)
            .add("mode_impute_acc", mode_impute_acc);
    }

    // Unmask x_val for eval.
    for (SparseX& s : x_val_true) {
        x_val[idx2d(s.i, s.l, HP.L)] = s.x;
    }

    // Tree eval.
    if (tree_fname != nullptr) {
        if ((variant_pos_fname == nullptr) || (variant_start_pos < 0)) {
            throw std::invalid_argument("Tree eval requires variant position file and variant start pos.");
        }
        std::vector<int> variant_pos{parse_pos_file_idx(variant_pos_fname, variant_start_pos, HP.L)};
        auto [coal_trees, recomb_pos] = parse_tree_file(tree_fname);
        std::vector<int> tree_idxs{get_tree_idxs(variant_pos, recomb_pos)};

        int excess_parsimony = 0;
        int emission_excess_parsimony = 0;
        std::chrono::steady_clock::duration t_parsimony_duration{};

        double weighted_clade_iou_sum = 0.0;
        double clade_weight_sum = 0.0;
        double emission_weighted_clade_iou_sum = 0.0;
        double emission_clade_weight_sum = 0.0;
        std::chrono::steady_clock::duration t_clade_iou_duration{};

        for (int l = 0; l < HP.L; ++l) {
            // Parismony.
            auto t_parsimony0 = std::chrono::steady_clock::now();

            std::unordered_map<Cluster*, int> cluster_idxs;
            cluster_idxs.reserve(clusters.rs[l].size());
            for (Cluster* c : clusters.rs[l]) {
                cluster_idxs.emplace(c, cluster_idxs.size());
            }
            std::vector<int> cluster_assign(HP.N);
            for (int i = 0; i < HP.N; ++i) {
                int idx = i < n_train_seqs ? train_idxs[i] : val_idxs[i-n_train_seqs];
                cluster_assign[idx] = cluster_idxs.at(clusters.r_assign[idx2d(i, l, HP.L)]);
            }
            excess_parsimony += calc_excess_parsimony(
                coal_trees[tree_idxs[l]], cluster_assign, clusters.rs[l].size()
            );

            std::vector<int> emission_counts(HP.K, 0);
            std::vector<int> emission_clusters(HP.N);
            for (int i = 0; i < HP.N; ++i) {
                int idx = i < n_train_seqs ? train_idxs[i] : val_idxs[i-n_train_seqs];
                int xi = comb_xi(x_train, x_val, idx2d(i,l,HP.L), n_train_seqs*HP.L);
                emission_clusters[idx] = xi;
                ++emission_counts[xi];
            }
            int n_obs_emissions = count_observed_labels(emission_clusters, HP.K);
            emission_excess_parsimony += calc_excess_parsimony(
                coal_trees[tree_idxs[l]], emission_clusters, n_obs_emissions 
            );

            t_parsimony_duration += std::chrono::steady_clock::now() - t_parsimony0;

            // Importance-weighted clade iou.
            auto t_clade_iou0 = std::chrono::steady_clock::now();
            for (const auto& [c, cluster_idx] : cluster_idxs) {
                double max_clade_iou = calc_max_clade_iou(
                    coal_trees[tree_idxs[l]], cluster_assign, cluster_idx, c->n
                );

                double z = static_cast<double>(c->n - 1) / (HP.N - 1);
                double weight = std::pow(z * (1.0 - z), clade_beta - 1.0);
                weighted_clade_iou_sum += weight * max_clade_iou;
                clade_weight_sum += weight;
            }

            for (int k = 0; k < HP.K; ++k) {
                if (emission_counts[k] == 0) {
                    continue;
                }
                double max_clade_iou = calc_max_clade_iou(
                    coal_trees[tree_idxs[l]], emission_clusters, k, emission_counts[k]
                );

                double z = static_cast<double>(emission_counts[k] - 1) / (HP.N - 1);
                double weight = std::pow(z * (1.0 - z), clade_beta - 1.0);
                emission_weighted_clade_iou_sum += weight * max_clade_iou;
                emission_clade_weight_sum += weight;
            }
            t_clade_iou_duration += std::chrono::steady_clock::now() - t_clade_iou0;

            // Tree viz.
            if ((tree_vis_fname != nullptr) && (l < 16)) {
                tree_to_dot(
                    tree_vis_fname, coal_trees[tree_idxs[l]], cluster_assign, emission_clusters, l, 16
                );
            }
        }
        double mean_excess_parsimony = static_cast<double>(excess_parsimony) / HP.L;
        double mean_emission_excess_parsimony = static_cast<double>(emission_excess_parsimony) / HP.L;

        double clade_iou = clade_weight_sum == 0.0 ? -1.0 : weighted_clade_iou_sum / clade_weight_sum;
        if (clade_weight_sum == 0.0) { std::cerr << "clade_iou undefined, no nontrivial clusters.\n"; }
        double emission_clade_iou = emission_clade_weight_sum == 0.0 ?
            -1.0 : emission_weighted_clade_iou_sum / emission_clade_weight_sum;
        if (emission_clade_weight_sum == 0.0) {std::cerr << "emission_clade_iou undefined, no nontrivial clusters.\n";}

        auto t_parsimony = std::chrono::duration_cast<std::chrono::milliseconds>(t_parsimony_duration).count();
        auto t_clade_iou = std::chrono::duration_cast<std::chrono::milliseconds>(t_clade_iou_duration).count();
        std::cerr << "mean_excess_parsimony=" << mean_excess_parsimony
            << " mean_emission_excess_parsimony=" << mean_emission_excess_parsimony
            << " t_parsimony=" << t_parsimony << "ms\n"
            << "clade_iou=" << clade_iou << " emission_clade_iou=" << emission_clade_iou
            << " clade_beta=" << clade_beta << " t_clade_iou=" << t_clade_iou << "ms\n";
        json.add("mean_excess_parsimony", mean_excess_parsimony)
            .add("mean_emission_excess_parsimony", mean_emission_excess_parsimony)
            .add("t_parsimony", t_parsimony)
            .add("clade_iou", clade_iou).add("emission_clade_iou", emission_clade_iou)
            .add("clade_beta", clade_beta).add("t_clade_iou", t_clade_iou);
    }

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

                int l_emission_same = comb_xi(x_train, x_val, idx2d(i,l,HP.L), n_train_seqs*HP.L)
                    == comb_xi(x_train, x_val, idx2d(j,l,HP.L), n_train_seqs*HP.L);
                int l1_emission_same = comb_xi(x_train, x_val, idx2d(i,l+1,HP.L), n_train_seqs*HP.L)
                    == comb_xi(x_train, x_val, idx2d(j,l+1,HP.L), n_train_seqs*HP.L);;
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
        cluster_purity /= HP.N*HP.L - n_masked_alleles;
    }
    else if (clusters.noisy) {
        cluster_purity = clusters.n_matches / static_cast<double>(clusters.n_obs);
    }
    std::cerr << "cluster_purity=" << cluster_purity << '\n';
    json.add("cluster_purity", cluster_purity);
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
    double val = -1.0;
    double mask = -1.0;

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

        else if (arg == "--val") { val = parse_double(argv[i+1]); }
        else if (arg == "--mask") { mask = parse_double(argv[i+1]); }

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
    bool do_val = val > 0.0;
    if (do_val != (mask > 0.0)) {
        throw std::invalid_argument("If imputation val frac > 0, need mask frac > 0.");
    };
    int n_train_seqs = 0;
    std::vector<int8_t> x_train;
    std::vector<size_t> train_idxs;
    x_train.reserve(static_cast<size_t>(HP.N*HP.L * (1.0-val)));
    train_idxs.reserve(static_cast<size_t>(HP.N*HP.L * (1.0-val)));

    int n_val_seqs = 0;
    std::vector<int8_t> x_val;
    std::vector<size_t> val_idxs;
    x_val.reserve(static_cast<size_t>(HP.N*HP.L * val));
    val_idxs.reserve(static_cast<size_t>(HP.N*HP.L * val));

    int n_masked_alleles = 0;
    std::vector<SparseX> x_val_true;
    x_val_true.reserve(static_cast<size_t>(HP.N*HP.L * val*mask));

    if (val > 0.0) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::bernoulli_distribution val_dist(val);
        std::bernoulli_distribution mask_dist(mask);

        for (int i = 0; i < HP.N; ++i) {
            auto line = x_raw.begin() + i*HP.L;
            if (!val_dist(gen)) {
                x_train.insert(x_train.end(), line, line+HP.L);
                train_idxs[n_train_seqs] = i;
                ++n_train_seqs;
                continue;
            }
            x_val.insert(x_val.end(), line, line+HP.L);
            val_idxs[n_val_seqs] = i;
            for (int l = 0; l < HP.L; ++l) {
                if (!mask_dist(gen)) { continue; }
                x_val_true.emplace_back(n_val_seqs, l, line[l]);
                x_val[idx2d(n_val_seqs, l, HP.L)] = -1;
                ++n_masked_alleles;
            }
            ++n_val_seqs;
        }
        HP.N = n_train_seqs;
        if (n_train_seqs == 0) { throw std::runtime_error("no train seqs."); }
        if (n_masked_alleles == 0) { throw std::runtime_error("invalid validation split."); }
        std::cerr << HP << "\nn_val_seqs=" << n_val_seqs << " n_masked_alleles=" << n_masked_alleles << '\n';
        json.add("n_val_seqs", n_val_seqs).add("n_masked_alleles", n_masked_alleles);
    }
    else {
        x_train = std::move(x_raw);
        for (int i = 0; i < HP.N; ++i) {
            train_idxs[i] = i;
        }
        n_train_seqs = HP.N;
    }

    dfcp(
        HP,

        soft, noisy,

        block_init,
        pbwt_init, pbwt_match_len, pbwt_match_curr,
        init_only,

        x_train, n_train_seqs, train_idxs,
        x_val, n_val_seqs, val_idxs,
        x_val_true, n_masked_alleles,

        tree_fname, variant_pos_fname, variant_start_pos,
        tree_vis_fname, clade_beta,
        json
    );
    std::cout << json.str() << '\n';
}

