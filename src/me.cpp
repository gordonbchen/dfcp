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
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "max.hpp"
#include "expect.hpp"
#include "elbo.hpp"
#include "tree.hpp"
#include "json.hpp"
#include "util.hpp"


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

    bool soft = false;
    bool block_init = false;

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

        else if (arg == "--val") { val = parse_double(argv[i+1]); }
        else if (arg == "--mask") { mask = parse_double(argv[i+1]); }

        else if (arg == "--tree") { tree_fname = argv[i+1]; }
        else if (arg == "--variant_pos") { variant_pos_fname = argv[i+1]; }
        else if (arg == "--variant_start_pos") { variant_start_pos = parse_int(argv[i+1]); }
        else if (arg == "--tree_vis") { tree_vis_fname = argv[i+1]; }

        else if (arg == "--soft") { soft = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--block_init") { block_init = (parse_int(argv[i+1]) == 1); }

        else { throw std::invalid_argument("Arg not recognized."); }
        i += 2;
    }
    std::cerr << HP << '\n';

    // Split val for imputation.
    bool do_val = val > 0.0;
    if (do_val != (mask > 0.0)) {
        throw std::invalid_argument("If imputation val frac > 0, need mask frac > 0.");
    };
    std::vector<int8_t> x(HP.N * HP.L, -1);
    int n_train_seqs = 0;
    int n_val_seqs = 0;
    std::vector<int> raw_to_split_idxs(HP.N, -1);

    int n_masked_alleles = 0;
    std::vector<SparseX> x_val_true;
    x_val_true.reserve(do_val ? static_cast<size_t>(val * mask * HP.N * HP.L) : 0);

    if (val > 0.0) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::bernoulli_distribution val_dist(val);
        std::bernoulli_distribution mask_dist(mask);

        for (int i = 0; i < HP.N; ++i) {
            auto line = x_raw.begin() + i*HP.L;
            if (!val_dist(gen)) {
                std::copy(line, line+HP.L, x.begin() + n_train_seqs*HP.L);
                raw_to_split_idxs[i] = n_train_seqs;
                ++n_train_seqs;
                continue;
            }
            std::copy(line, line+HP.L, x.end() - (n_val_seqs+1)*HP.L);
            raw_to_split_idxs[i] = HP.N - (n_val_seqs+1);
            for (int l = 0; l < HP.L; ++l) {
                if (!mask_dist(gen)) { continue; }
                x_val_true.emplace_back(HP.N - (n_val_seqs+1), l, line[l]);
                x[idx2d(HP.N - (n_val_seqs+1), l, HP.L)] = -1;
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
        x = std::move(x_raw);
        for (int i = 0; i < HP.N; ++i) {
            raw_to_split_idxs[i] = i;
        }
        n_train_seqs = HP.N;
    }

    // Init params and clusters.
    Params params{HP};
    Clusters clusters{HP, soft};
    if (block_init) {
        clusters.block_init(x);
    }
    else {
        HP.N = 0;
        add_seqs(clusters, x.begin(), n_train_seqs, HP, params);
    }

    EarlyStopping early_stop{2, false, 1e-3};
    double elbo = 0.0;

    std::vector<Json> train_log;
    while (!early_stop.converged()) {
        auto t0 = std::chrono::steady_clock::now();
        max_step(clusters, x, HP, params);
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
    param_log.add("mu_alpha", params.mu_alpha);
    param_log.add("mu_gamma", params.mu_gamma);
    param_log.add("mu_d", params.mu_d);
    json.add("params", param_log);

    // Impute.
    if (n_val_seqs > 0) {
        auto t0 = std::chrono::steady_clock::now();
        add_seqs(clusters, x.begin() + n_train_seqs*HP.L, n_val_seqs, HP, params);
        auto t1 = std::chrono::steady_clock::now();
        auto t_impute = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        std::vector<int8_t> modes{count_modes(x, n_train_seqs, HP.L, HP.K)};

        int n_dfcp_correct = 0;
        int n_mode_correct = 0;
        for (SparseX& s : x_val_true) {
            if (s.x == clusters.r_assign[idx2d(s.i, s.l, HP.L)]->get_imputed_emission()) {
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

    // Unmask x for eval.
    for (SparseX& s : x_val_true) {
        x[idx2d(s.i, s.l, HP.L)] = s.x;
    }

    // Tree parsimony.
    if (tree_fname != nullptr) {
        if ((variant_pos_fname == nullptr) || (variant_start_pos < 0)) {
            throw std::invalid_argument("Tree eval requires variant position file and variant start pos.");
        }
        std::vector<int> variant_pos{parse_pos_file_idx(variant_pos_fname, variant_start_pos, HP.L)};
        auto [coal_trees, recomb_pos] = parse_tree_file(tree_fname);
        std::vector<int> tree_idxs{get_tree_idxs(variant_pos, recomb_pos)};

        auto t0 = std::chrono::steady_clock::now();
        int excess_parsimony = 0;
        int emission_excess_parsimony = 0;
        for (int l = 0; l < HP.L; ++l) {
            std::unordered_map<Cluster*, int> cluster_idxs;
            cluster_idxs.reserve(clusters.rs[l].size());
            for (Cluster* c : clusters.rs[l]) {
                cluster_idxs.emplace(c, cluster_idxs.size());
            }
            std::vector<int> cluster_assign(HP.N);
            for (int i = 0; i < HP.N; ++i) {
                cluster_assign[i] = cluster_idxs.at(clusters.r_assign[idx2d(raw_to_split_idxs[i], l, HP.L)]);
            }
            excess_parsimony += calc_excess_parsimony(
                coal_trees[tree_idxs[l]], cluster_assign, clusters.rs[l].size()
            );

            std::vector<int> emission_clusters(HP.N);
            for (int i = 0; i < HP.N; ++i) {
                emission_clusters[i] = x[idx2d(raw_to_split_idxs[i], l, HP.L)];
            }
            int n_obs_emissions = count_observed_labels(emission_clusters, HP.K);
            emission_excess_parsimony += calc_excess_parsimony(
                coal_trees[tree_idxs[l]], emission_clusters, n_obs_emissions 
            );

            if ((tree_vis_fname != nullptr) && (l < 8)) {
                tree_to_dot(
                    tree_vis_fname, coal_trees[tree_idxs[l]], cluster_assign, emission_clusters, l, 8
                );
            }
        }
        double mean_excess_parsimony = static_cast<double>(excess_parsimony) / HP.L;
        double mean_emission_excess_parsimony = static_cast<double>(emission_excess_parsimony) / HP.L;

        auto t1 = std::chrono::steady_clock::now();
        auto t_parsimony = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        std::cerr << "mean_excess_parsimony=" << mean_excess_parsimony
            << " mean_emission_excess_parsimony=" << mean_emission_excess_parsimony
            << " t_parsimony=" << t_parsimony << "ms\n";
        json.add("mean_excess_parsimony", mean_excess_parsimony)
            .add("mean_emission_excess_parsimony", mean_emission_excess_parsimony)
            .add("t_parsimony", t_parsimony);
    }

    // IOU.
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

                int l_emission_same = x[idx2d(i,l,HP.L)] == x[idx2d(j,l,HP.L)];
                int l1_emission_same = x[idx2d(i,l+1,HP.L)] == x[idx2d(j,l+1,HP.L)];
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

    std::cout << json.str() << '\n';
    return 0;
}

