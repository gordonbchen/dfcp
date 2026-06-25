#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <vector>
#include <cstring>
#include "tree.hpp"
#include "clusters.hpp"
#include "util.hpp"


std::vector<int> parse_pos_file_idx(const char *fname, int start_idx, int L) {
    std::ifstream pos_file{fname};
    if (!pos_file.is_open()) { throw std::invalid_argument("Failed to open position file."); }

    int n_pos;
    pos_file >> n_pos;

    if (start_idx + L > n_pos) { throw std::invalid_argument("Position file not long enough for start and length"); };

    std::vector<int> pos;
    pos.reserve(L);

    int x;
    for (int i = 0; i < start_idx; ++i) {
        pos_file >> x;
        if (pos_file.peek() == ',') { pos_file.ignore(); }
    }

    for (int i = 0; i < L; ++i) {
        pos_file >> x;
        pos.push_back(x);
        if (pos_file.peek() == ',') { pos_file.ignore(); }
    }
    return pos;
}

const char* parse_coal_subtree(const char *s, std::unordered_map<int, std::tuple<int, int>>& coal_tree, int idx) {
    while (*s != '(') { ++s; }
    ++s;

    int left_idx;
    if (*s == '(') {
        left_idx = -(2 * std::abs(idx));
        s = parse_coal_subtree(s, coal_tree, left_idx);
    }
    else {
        left_idx = std::strtol(s, nullptr, 10) - 1;
    }

    while (*s != ' ') { ++s; }
    ++s;

    int right_idx;
    if (*s == '(') {
        right_idx = -(2*std::abs(idx) + 1);
        s = parse_coal_subtree(s, coal_tree, right_idx);
    }
    else {
        right_idx = std::strtol(s, nullptr, 10) - 1;
    }

    coal_tree.emplace(idx, std::tuple<int, int>{left_idx, right_idx});
    return s;
}

int parse_coal_tree(const char *s, std::unordered_map<int, std::tuple<int, int>>& coal_tree) {
    s = std::strstr(s, "pos_");
    if (s == nullptr) { throw std::runtime_error("Failed to parse tree, could not find 'pos_'"); }
    s += 4;

    int recomb_pos = std::strtol(s, nullptr, 10);
    parse_coal_subtree(s, coal_tree, -1);
    return recomb_pos;
}

std::tuple<std::vector<std::unordered_map<int, std::tuple<int, int>>>, std::vector<int>> parse_tree_file(
    const char *fname
) {
    std::ifstream tree_file{fname};
    if (!tree_file.is_open()) { throw std::runtime_error("Failed to open tree file."); };

    std::string line;
    for (int i = 0; i < 3; ++i) {
        std::getline(tree_file, line);
    }

    std::vector<std::unordered_map<int, std::tuple<int, int>>> trees;
    std::vector<int> recomb_pos;
    int l = 0;
    while (std::getline(tree_file, line)) {
        if (line[0] == 'e') { break; }
        trees.emplace_back();
        recomb_pos.push_back(parse_coal_tree(line.c_str(), trees[l]));
        ++l;
    }
    return std::tuple{trees, recomb_pos};
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

