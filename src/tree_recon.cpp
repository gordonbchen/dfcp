#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <format>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <utility>
#include <CLI/CLI.hpp>
#include "r_assign_io.hpp"
#include "tree.hpp"
#include "util.hpp"

namespace {

struct Node {
    const BitSet* left;
    const BitSet* right;
};

std::pair<const BitSet*, const BitSet*> get_ord_pair(const BitSet *a, const BitSet *b) {
    return std::pair<const BitSet*,const BitSet*>{std::min(a, b), std::max(a, b)};
}

void to_coal_tree_recur(
    int curr_int, const BitSet* curr_bs,
    std::unordered_map<int, CoalNode>& coal_tree,
    int& next_internal_id,
    std::unordered_map<BitSet, Node, BitSetHash>& tree,
    std::unordered_map<const BitSet*, int>& leaves_map
) {
    Node& node = tree[*curr_bs];
    int left_int = leaves_map.contains(node.left) ? leaves_map[node.left] : next_internal_id--;
    int right_int = leaves_map.contains(node.right) ? leaves_map[node.right] : next_internal_id--;
    coal_tree[curr_int] = CoalNode{right_int, left_int, 0.0, 0.0};

    if (left_int < 0) {
        to_coal_tree_recur(left_int, node.left, coal_tree, next_internal_id, tree, leaves_map);
    }
    if (right_int < 0) {
        to_coal_tree_recur(right_int, node.right, coal_tree, next_internal_id, tree, leaves_map);
    }
}

std::unordered_map<int, CoalNode> to_coal_tree(
    const BitSet* root,
    std::unordered_map<BitSet, Node, BitSetHash>& tree,
    std::unordered_map<const BitSet*, int>& leaves_map
) {
    std::unordered_map<int, CoalNode> coal_tree;
    int next_internal_id = -2;
    to_coal_tree_recur(-1, root, coal_tree, next_internal_id, tree, leaves_map);
    return coal_tree;
}

void upgma(
    std::vector<RAssign>& r_assigns, int N, int l,
    int n_tree_viz, const std::optional<std::string>& recon_tree_viz
) {
    // Populate the tree with individual leaves.
    std::unordered_map<BitSet, Node, BitSetHash> tree;
    std::vector<const BitSet*> leaves;
    std::unordered_map<const BitSet*, int> leaves_map;
    std::unordered_set<const BitSet*> active_nodes;
    for (int i = 0; i < N; ++i) {
        BitSet bs{static_cast<size_t>(N)};
        bs.set(i);
        const BitSet *bs_ptr = &(tree.emplace(std::move(bs), Node{}).first->first);
        leaves.push_back(bs_ptr);
        leaves_map[bs_ptr] = i;
        active_nodes.emplace(bs_ptr);
    }

    // Compute pairwise scores.
    std::unordered_map<std::pair<const BitSet*, const BitSet*>, double, PairPointerHash> scores;
    for (RAssign& r_assign : r_assigns) {
        std::unordered_map<uint32_t, std::unordered_set<int>> clusters;
        for (int i = 0; i < N; ++i) {
            auto it = clusters.find(r_assign(i,l));
            if (it == clusters.end()) {
                it = clusters.emplace(r_assign(i,l), std::unordered_set<int>()).first;
            }

            for (int j = 0; j < i; ++j) {
                scores[get_ord_pair(leaves[i], leaves[j])] += it->second.contains(j);
            }
            it->second.emplace(i);
        }
    }

    while (active_nodes.size() > 1) {
        // Find the nodes with the highest pair score;
        auto best_pair = std::max_element(
            scores.begin(), scores.end(),
            [](const auto& p1, const auto& p2) { return p1.second < p2.second; }
        );
        const BitSet *left_bs = best_pair->first.first;
        const BitSet *right_bs = best_pair->first.second;

        // Add new node to the tree.
        BitSet new_bs{static_cast<size_t>(N)};
        new_bs.sunion(*left_bs, *right_bs);
        const BitSet *new_bs_ptr = &(tree.emplace(std::move(new_bs), Node{left_bs, right_bs}).first->first);

        // Remove joined nodes.
        active_nodes.erase(left_bs);
        active_nodes.erase(right_bs);
        scores.erase(best_pair);

        // Update old node scores.
        for (const BitSet* node : active_nodes) {
            auto left_it = scores.find(get_ord_pair(left_bs, node));
            auto right_it = scores.find(get_ord_pair(right_bs, node));
            if ((left_it == scores.end()) || (right_it == scores.end())) {
                throw std::runtime_error("Could not find removed bitset and leaf in scores.");
            }
            double new_score = (
                ((left_bs->n_set * left_it->second) + (right_bs->n_set * right_it->second)) /
                (left_bs->n_set + right_bs->n_set)
            );
            scores.erase(left_it);
            scores.erase(right_it);
            scores[get_ord_pair(new_bs_ptr, node)] = new_score;
        }

        // Add new node.
        active_nodes.emplace(new_bs_ptr);
    }

    // Tree viz.
    if (recon_tree_viz && (l < n_tree_viz)) {
        const BitSet *root = *active_nodes.begin();
        std::unordered_map<int, CoalNode> coal_tree{to_coal_tree(root, tree, leaves_map)};
        tree_to_dot(*recon_tree_viz, coal_tree, std::nullopt, std::nullopt, l, n_tree_viz);
    }
}
}

int main(int argc, char *argv[]) {
    // Parse args.
    CLI::App app{"Tree Recon"};

    std::vector<std::string> cluster_files;
    app.add_option("cluster_files", cluster_files, "Cluster assignment bins.")
        ->expected(1, -1)->required();

    std::optional<std::string> true_trees;
    app.add_option("--true_trees", true_trees, "True tree files.");

    std::optional<std::string> recon_tree_viz;
    app.add_option("--recon_tree_viz", recon_tree_viz, "Recon tree viz file.");

    CLI11_PARSE(app, argc, argv);

    // Read r_assigns.
    std::vector<RAssign> r_assigns;
    int N = -1;
    int L = -1;
    for (std::string& file : cluster_files) {
        std::cerr << "Reading " << file << ".\n";
        r_assigns.push_back(read_r_assign(file));
        if (N == -1) {
            N = r_assigns.back().N;
            L = r_assigns.back().L;
            continue;
        }
        if ((r_assigns.back().N != N) || (r_assigns.back().L != L)) {
            throw std::runtime_error(
                std::format("Expected N={} L={}, got N={} L={} from file '{}'",
                            N, L, r_assigns.back().N, r_assigns.back().L, file)
            );
        }
    }

    std::cerr << "Running UPGMA.\n";
    for (int l = 0; l < L; ++l) {
        upgma(r_assigns, N, l, std::min(16, L), recon_tree_viz);
    }
}
