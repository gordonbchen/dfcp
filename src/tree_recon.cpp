#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <format>
#include <stdexcept>
#include <unordered_map>
#include <optional>
#include <CLI/CLI.hpp>
#include "r_assign_io.hpp"
#include "tree.hpp"

namespace {

uint64_t get_ord_pair(int a, int b) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(std::min(a, b))) << 32)
         | static_cast<uint32_t>(std::max(a, b));
}

std::unordered_map<int, CoalNode> upgma(const std::vector<RAssign>& r_assigns, int N, int l) {
    std::unordered_map<int, CoalNode> tree;
    std::unordered_map<int, int> sizes;
    for (int i = 0; i < N; ++i) {
        sizes.emplace(i, 1);
    }

    std::unordered_map<uint64_t, double> scores;
    for (const RAssign& r_assign : r_assigns) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < i; ++j) {
                scores[get_ord_pair(i, j)] += r_assign(i, l) == r_assign(j, l);
            }
        }
    }

    // Allocate internal IDs in merge order so the final root is -1.
    int next_internal_id = -(N - 1);
    while (sizes.size() > 1) {
        auto best_pair = std::max_element(
            scores.begin(), scores.end(),
            [](const auto& p1, const auto& p2) { return p1.second < p2.second; }
        );
        int left = static_cast<int32_t>(best_pair->first >> 32);
        int right = static_cast<int32_t>(static_cast<uint32_t>(best_pair->first));
        int left_size = sizes.at(left);
        int right_size = sizes.at(right);
        int node = next_internal_id++;
        tree.emplace(node, CoalNode{left, right, 0.0, 0.0});

        sizes.erase(left);
        sizes.erase(right);
        scores.erase(best_pair);
        for (const auto& [other, size] : sizes) {
            auto left_it = scores.find(get_ord_pair(left, other));
            auto right_it = scores.find(get_ord_pair(right, other));
            if (left_it == scores.end() || right_it == scores.end()) {
                throw std::runtime_error("Could not find merged node pair in scores.");
            }
            double score = (left_size * left_it->second + right_size * right_it->second)
                         / (left_size + right_size);
            scores.erase(left_it);
            scores.erase(right_it);
            scores.emplace(get_ord_pair(node, other), score);
        }
        sizes.emplace(node, left_size + right_size);
    }
    return tree;
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
        auto tree = upgma(r_assigns, N, l);
        int n_tree_viz = std::min(16, L);
        if (recon_tree_viz && l < n_tree_viz) {
            tree_to_dot(*recon_tree_viz, tree, std::nullopt, std::nullopt, l, n_tree_viz);
        }
    }
}
