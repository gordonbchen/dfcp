#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <string_view>
#include <cstdlib>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "tree.hpp"
#include "json.hpp"
#include "util.hpp"


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


std::pair<std::vector<int8_t>, std::vector<int>> read_seq_file(char *filename, int& N, int *L, bool allow_missing) {
    std::ifstream file(filename);
    if (!file.is_open()) { throw std::runtime_error("Failed to open seq file."); };

    std::vector<int8_t> x;
    std::vector<int> masked_ls;

    std::string line;
    N = 0;
    while (std::getline(file, line)) {
        if ((L == nullptr) && (N == 0)) {
            *L = line.size();
        }

        for (size_t l = 0; l < line.size(); ++l) {
            if (allow_missing && (line[l] == '.')) {
                x.push_back(-1);
                if (N == 0) {
                    masked_ls.push_back(l);
                }
            }
            else if ((line[l] >= '0') && (line[l] <= '9')) {
                x.push_back(line[l] - '0');
            }
            else {
                throw std::runtime_error("Invalid allele char.");
            };
        }
        ++N;
    }
    return {std::move(x), std::move(masked_ls)};
}

int main(int argc, char *argv[]) {
    Json json;

    // Read ref and target files.
    if (argc < 3) { throw std::invalid_argument("Requires reference and target seq files."); }
    std::cerr << "ref_file=" << argv[1] << " target_file=" << argv[2] << '\n';
    json.add("ref_file", argv[1]).add("target_file", argv[2]);

    int n_train_seqs;
    int L;
    std::vector<int8_t> x_train{read_seq_file(argv[1], n_train_seqs, &L, false).first};

    int n_val_seqs;
    auto [x_val, masked_ls] = read_seq_file(argv[2], n_val_seqs, nullptr, true);

    HyperParams HP{.N=n_train_seqs, .L=L};

    // Parse optional args.
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
            x_val_true, masked_ls,

            json
        );
    }
    std::cout << json.str() << '\n';
}
