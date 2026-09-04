#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include "json.hpp"
#include "r_assign_io.hpp"
#include "seq_array.hpp"
#include "tree.hpp"
#include "util.hpp"


namespace {

std::uint64_t choose_two(std::uint64_t n) {
    return n * (n - 1) / 2;
}

std::pair<double, double> adjacent_ious(const RAssign& r_assign, const SeqArray& x) {
    double cluster_sum = 0.0;
    double emission_sum = 0.0;
    for (int l = 0; l < r_assign.L - 1; ++l) {
        std::unordered_map<std::uint32_t, std::uint64_t> left;
        std::unordered_map<std::uint32_t, std::uint64_t> right;
        std::unordered_map<std::uint64_t, std::uint64_t> joint;
        std::array<std::uint64_t, 2> emission_left{};
        std::array<std::uint64_t, 2> emission_right{};
        std::array<std::uint64_t, 4> emission_joint{};

        for (int i = 0; i < r_assign.N; ++i) {
            std::uint32_t a = r_assign(i, l);
            std::uint32_t b = r_assign(i, l + 1);
            ++left[a];
            ++right[b];
            ++joint[(static_cast<std::uint64_t>(a) << 32) | b];

            int left_x = x(i, l);
            int right_x = x(i, l + 1);
            ++emission_left[left_x];
            ++emission_right[right_x];
            ++emission_joint[2 * left_x + right_x];
        }

        std::uint64_t intersection = 0;
        std::uint64_t left_pairs = 0;
        std::uint64_t right_pairs = 0;
        for (const auto& [key, count] : joint) { intersection += choose_two(count); }
        for (const auto& [key, count] : left) { left_pairs += choose_two(count); }
        for (const auto& [key, count] : right) { right_pairs += choose_two(count); }
        std::uint64_t union_size = left_pairs + right_pairs - intersection;
        cluster_sum += union_size == 0 ? 1.0 : static_cast<double>(intersection) / union_size;

        std::uint64_t emission_intersection = 0;
        std::uint64_t emission_left_pairs = 0;
        std::uint64_t emission_right_pairs = 0;
        for (std::uint64_t count : emission_joint) { emission_intersection += choose_two(count); }
        for (std::uint64_t count : emission_left) { emission_left_pairs += choose_two(count); }
        for (std::uint64_t count : emission_right) { emission_right_pairs += choose_two(count); }
        std::uint64_t emission_union = emission_left_pairs + emission_right_pairs
            - emission_intersection;
        emission_sum += emission_union == 0 ? 1.0
            : static_cast<double>(emission_intersection) / emission_union;
    }
    return {cluster_sum / (r_assign.L - 1), emission_sum / (r_assign.L - 1)};
}

void eval_partitions(const RAssign& r_assign, const SeqArray& x, Json& json) {
    auto [mean_adj_iou, mean_adj_emission_iou] = adjacent_ious(r_assign, x);

    std::uint64_t n_clusters = 0;
    std::uint64_t majority_alleles = 0;
    for (int l = 0; l < r_assign.L; ++l) {
        std::unordered_map<std::uint32_t, std::array<int, 2>> counts;
        for (int i = 0; i < r_assign.N; ++i) {
            ++counts[r_assign(i, l)][x(i, l)];
        }
        n_clusters += counts.size();
        for (const auto& [id, count] : counts) {
            majority_alleles += std::max(count[0], count[1]);
        }
    }

    double mean_clusters = static_cast<double>(n_clusters) / r_assign.L;
    double purity = static_cast<double>(majority_alleles) / r_assign.ids.size();
    std::cerr << "mean_adj_iou=" << mean_adj_iou << " mean_adj_emission_iou=" << mean_adj_emission_iou
        << " mean_clusters=" << mean_clusters << " cluster_purity=" << purity << '\n';
    json.add("mean_adj_iou", mean_adj_iou).add("mean_adj_emission_iou", mean_adj_emission_iou)
        .add("mean_clusters", mean_clusters).add("cluster_purity", purity);
}


double clade_weight(int cluster_size, int n_sequences, double beta) {
    double z = static_cast<double>(cluster_size - 1) / (n_sequences - 1);
    return std::pow(z * (1.0 - z), beta - 1.0);
}

void eval_trees(
    const RAssign& r_assign, const SeqArray& x,
    const std::vector<int>& variant_pos, const char* tree_file,
    double clade_beta, const char* tree_vis_file, Json& json
) {
    auto [trees, recomb_pos] = parse_tree_file(tree_file);
    std::vector<int> tree_idxs = get_tree_idxs(variant_pos, recomb_pos);

    std::vector<bool> tree_leaves(r_assign.N);
    for (const auto& [idx, node] : trees.front()) {
        for (int child : {node.left, node.right}) {
            if (child >= r_assign.N) {
                throw std::runtime_error("Tree leaf is outside the reference sequence array.");
            }
            if (child >= 0) { tree_leaves[child] = true; }
        }
    }
    if (std::find(tree_leaves.begin(), tree_leaves.end(), false) != tree_leaves.end()) {
        throw std::runtime_error("Tree does not contain every reference sequence.");
    }

    std::uint64_t excess_parsimony = 0;
    std::uint64_t emission_excess_parsimony = 0;
    double weighted_clade_iou = 0.0;
    double clade_weight_sum = 0.0;
    double weighted_emission_clade_iou = 0.0;
    double emission_clade_weight_sum = 0.0;

    std::vector<int> dense_r_assign(r_assign.N);
    std::vector<int> emissions(r_assign.N);
    int visualized_loci = std::min(r_assign.L, 16);
    for (int l = 0; l < r_assign.L; ++l) {
        std::unordered_map<std::uint32_t, int> dense_ids;
        std::vector<int> cluster_sizes;
        std::array<int, 2> emission_counts{};
        for (int i = 0; i < r_assign.N; ++i) {
            auto [it, inserted] = dense_ids.try_emplace(r_assign(i, l), dense_ids.size());
            if (inserted) { cluster_sizes.push_back(0); }
            dense_r_assign[i] = it->second;
            ++cluster_sizes[it->second];
            emissions[i] = x(i, l);
            ++emission_counts[emissions[i]];
        }

        const std::unordered_map<int, CoalNode>& tree = trees[tree_idxs[l]];
        int n_clusters = static_cast<int>(cluster_sizes.size());
        excess_parsimony += calc_excess_parsimony(
            tree, dense_r_assign, n_clusters, n_clusters
        );
        int n_emissions = (emission_counts[0] != 0) + (emission_counts[1] != 0);
        emission_excess_parsimony += calc_excess_parsimony(
            tree, emissions, n_emissions, 2
        );

        for (int cluster = 0; cluster < n_clusters; ++cluster) {
            double iou = calc_max_clade_iou(
                tree, dense_r_assign, cluster, cluster_sizes[cluster]
            ).first;
            double weight = clade_weight(cluster_sizes[cluster], r_assign.N, clade_beta);
            weighted_clade_iou += weight * iou;
            clade_weight_sum += weight;
        }
        for (int allele = 0; allele < 2; ++allele) {
            if (emission_counts[allele] == 0) { continue; }
            double iou = calc_max_clade_iou(
                tree, emissions, allele, emission_counts[allele]
            ).first;
            double weight = clade_weight(emission_counts[allele], r_assign.N, clade_beta);
            weighted_emission_clade_iou += weight * iou;
            emission_clade_weight_sum += weight;
        }

        if (tree_vis_file != nullptr && l < visualized_loci) {
            tree_to_dot(
                tree_vis_file, tree, dense_r_assign, emissions, l, visualized_loci
            );
        }
    }

    double mean_excess_parsimony = static_cast<double>(excess_parsimony) / r_assign.L;
    double mean_emission_excess_parsimony =
        static_cast<double>(emission_excess_parsimony) / r_assign.L;
    double clade_iou = clade_weight_sum == 0.0 ? -1.0 : weighted_clade_iou / clade_weight_sum;
    double emission_clade_iou = emission_clade_weight_sum == 0.0 ? -1.0
        : weighted_emission_clade_iou / emission_clade_weight_sum;

    std::cerr << "mean_excess_parsimony=" << mean_excess_parsimony
        << " mean_emission_excess_parsimony=" << mean_emission_excess_parsimony
        << " clade_iou=" << clade_iou << " emission_clade_iou=" << emission_clade_iou << '\n';
    json.add("mean_excess_parsimony", mean_excess_parsimony)
        .add("mean_emission_excess_parsimony", mean_emission_excess_parsimony)
        .add("clade_iou", clade_iou).add("emission_clade_iou", emission_clade_iou)
        .add("clade_beta", clade_beta);
}

}


int main(int argc, char* argv[]) {
    if (argc < 5) {
        throw std::invalid_argument(
            "Requires reference, R assignments, variant positions, and tree files."
        );
    }

    double clade_beta = 2.0;
    const char* tree_vis_file = nullptr;
    for (int i = 5; i < argc; i += 2) {
        if (i + 1 >= argc) { throw std::invalid_argument("Arg has no value."); }
        std::string_view arg{argv[i]};

        if (arg == "--clade_beta") { clade_beta = parse_double(argv[i + 1]); }
        else if (arg == "--tree_vis") { tree_vis_file = argv[i + 1]; }

        else { throw std::invalid_argument("Arg not recognized."); }
    }
    if (clade_beta < 1.0) { throw std::invalid_argument("clade_beta must be at least 1."); }

    SeqArray x = read_seq_file(argv[1]);
    RAssign r_assign = read_r_assign(argv[2]);
    std::vector<int> variant_pos = read_variant_pos(argv[3]);
    if (r_assign.N != x.N || r_assign.L != x.L || static_cast<int>(variant_pos.size()) != x.L) {
        throw std::runtime_error(
            "Reference, R assignments, and variant positions have different dimensions."
        );
    }

    Json json;
    json.add("ref_file", argv[1]).add("r_assign_file", argv[2])
        .add("variant_pos_file", argv[3]).add("tree_file", argv[4]);

    eval_partitions(r_assign, x, json);
    eval_trees(r_assign, x, variant_pos, argv[4], clade_beta, tree_vis_file, json);

    std::cout << json.str() << '\n';
}
