#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <vector>
#include "tree.hpp"
#include "clusters.hpp"
#include "util.hpp"


const char* parse_coal_tree(const char *s, std::unordered_map<int, std::tuple<int, int>>& coal_tree, int idx) {
    while (*s != '(') { ++s; }
    ++s;

    int left_idx;
    if (*s == '(') {
        left_idx = -(2 * std::abs(idx));
        s = parse_coal_tree(s, coal_tree, left_idx);
    }
    else {
        left_idx = std::strtol(s, nullptr, 10) - 1;
    }

    while (*s != ' ') { ++s; }
    ++s;

    int right_idx;
    if (*s == '(') {
        right_idx = -(2*std::abs(idx) + 1);
        s = parse_coal_tree(s, coal_tree, right_idx);
    }
    else {
        right_idx = std::strtol(s, nullptr, 10) - 1;
    }

    coal_tree.emplace(idx, std::tuple<int, int>{left_idx, right_idx});
    return s;
}

std::vector<std::unordered_map<int, std::tuple<int, int>>> parse_tree_file(const char *tree_file_name, int L) {
    std::ifstream tree_file{tree_file_name};
    if (!tree_file.is_open()) { throw std::runtime_error("Failed to open tree file."); };

    std::string line;
    for (int i = 0; i < 3; ++i) {
        std::getline(tree_file, line);
    }

    std::vector<std::unordered_map<int, std::tuple<int, int>>> trees(L);
    int l = 0;
    while (std::getline(tree_file, line)) {
        if (line[0] == 'e') { break; }
        if (l >= L) { throw std::runtime_error("# trees should match sequence length"); }
        parse_coal_tree(line.c_str(), trees[l], -1);
        ++l;
    }
    return trees;
}

struct ParsimonyMsg {
    int score;
    uint64_t cluster_bm;
};

ParsimonyMsg calc_parsimony(
    int idx,
    int l,
    const std::unordered_map<int, std::tuple<int, int>>& coal_tree,
    const Clusters& clusters,
    const std::unordered_map<Cluster*, int>& cluster_idxs
) {
    if (clusters.rs[l].size() > 64) {
        throw std::runtime_error("uint64_t insufficient for cluster set bitvector.");
    };
    uint64_t cluster_bm = 0;
    if (!coal_tree.contains(idx)) {
        cluster_bm |= 1ULL << cluster_idxs.at(clusters.r_assign[idx2d(idx, l, clusters.HP.L)]);
        return ParsimonyMsg{0, cluster_bm};
    }

    const auto& [left, right] = coal_tree.at(idx);
    ParsimonyMsg lp = calc_parsimony(left, l, coal_tree, clusters, cluster_idxs);
    ParsimonyMsg rp = calc_parsimony(right, l, coal_tree, clusters, cluster_idxs);

    cluster_bm = lp.cluster_bm & rp.cluster_bm;
    if (cluster_bm != 0) {
        return ParsimonyMsg{lp.score + rp.score, cluster_bm};
    }
    return ParsimonyMsg{lp.score + rp.score + 1, lp.cluster_bm | rp.cluster_bm};
}

int calc_excess_parsimony(
    int l,
    const std::unordered_map<int, std::tuple<int, int>>& coal_tree,
    const Clusters& clusters,
    const std::unordered_map<Cluster*, int>& cluster_idxs
) {
    int parsimony = calc_parsimony(-1, l, coal_tree, clusters, cluster_idxs).score;
    return parsimony - (clusters.rs[l].size() - 1);
}

