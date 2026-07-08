#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ios>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <utility>
#include <vector>
#include <cstring>
#include <bit>
#include <format>
#include "tree.hpp"


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


struct BitSet {
    size_t size;
    size_t n_words;
    std::vector<uint64_t> bm;

    BitSet(size_t size_) : size(size_), n_words((size_ + 63) >> 6), bm((size_ + 63) >> 6, 0) {}

    void set(size_t idx) {
        if (idx >= size) {
            throw std::invalid_argument("idx is too large for bit set.");
        }
        bm[idx >> 6] |= 1ULL << (((1 << 6) - 1) & idx);
    }

    bool is_empty() {
        for (size_t i = 0; i < n_words; ++i) {
            if (bm[i] != 0) {
                return false;
            }
        }
        return true;
    }

    void intersect(const BitSet& l, const BitSet& r) {
        if (!((size == l.size) && (l.size == r.size))) {
            throw std::invalid_argument("BitSets must be the same size.");
        }
        for (size_t i = 0; i < n_words; ++i) {
            bm[i] = l.bm[i] & r.bm[i];
        }
    }

    void sunion(const BitSet& l, const BitSet& r) {
        if (!((size == l.size) && (l.size == r.size))) {
            throw std::invalid_argument("BitSets must be the same size.");
        }
        for (size_t i = 0; i < n_words; ++i) {
            bm[i] = l.bm[i] | r.bm[i];
        }
    }

    size_t count1s() {
        size_t n = 0;
        for (size_t i = 0; i < n_words; ++i) {
            n += std::popcount(bm[i]);
        }
        return n;
    }
};

int count_observed_labels(const std::vector<int>& labels, size_t max_size) {
    BitSet bm{max_size};
    for (int x : labels) {
        if (x < 0 || static_cast<size_t>(x) >= max_size) {
            throw std::runtime_error("label not b/t 0, max_size.");
        }
        bm.set(x);
    }
    return bm.count1s();
}

struct ParsimonyMsg {
    int score;
    BitSet cluster_bm;
};

ParsimonyMsg calc_parsimony(
    int idx,
    const std::unordered_map<int, std::tuple<int, int>>& coal_tree,
    const std::vector<int>& cluster_assignments,
    const size_t n_clusters
) {
    BitSet cluster_bm{n_clusters};
    if (!coal_tree.contains(idx)) {
        cluster_bm.set(cluster_assignments[idx]);
        return ParsimonyMsg{0, std::move(cluster_bm)};
    }

    const auto& [left, right] = coal_tree.at(idx);
    ParsimonyMsg lp = calc_parsimony(left, coal_tree, cluster_assignments, n_clusters);
    ParsimonyMsg rp = calc_parsimony(right, coal_tree, cluster_assignments, n_clusters);

    cluster_bm.intersect(lp.cluster_bm, rp.cluster_bm);
    if (!cluster_bm.is_empty()) {
        return ParsimonyMsg{lp.score + rp.score, std::move(cluster_bm)};
    }
    cluster_bm.sunion(lp.cluster_bm, rp.cluster_bm);
    return ParsimonyMsg{lp.score + rp.score + 1, std::move(cluster_bm)};
}

int calc_excess_parsimony(
    const std::unordered_map<int, std::tuple<int, int>>& coal_tree,
    const std::vector<int>& cluster_assignments,
    int n_clusters
) {
    int parsimony = calc_parsimony(-1, coal_tree, cluster_assignments, n_clusters).score;
    int excess_parsimony = parsimony - (n_clusters - 1);
    if (excess_parsimony < 0) { throw std::runtime_error("Negative excess parsimony."); }
    return excess_parsimony;
}

std::vector<int> get_tree_idxs(const std::vector<int>& variant_pos, const std::vector<int>& recomb_pos) {
    int L = variant_pos.size();
    std::vector<int> tree_idxs(L);
    int tree_idx = 0;
    for (int l = 0; l < L; ++l) {
        while ((tree_idx < static_cast<int>(recomb_pos.size()) - 1) && (recomb_pos[tree_idx+1] <= variant_pos[l])) {
            ++tree_idx;
        }
        tree_idxs[l] = tree_idx;
    }
    return tree_idxs;
}


std::string node_name(int x, int l) {
    return std::format("\"l{}_{}\"", l, x);
}

const std::vector<std::string>colors = {"yellowgreen", "cyan", "orange", "deeppink", "red"};
const std::vector<std::string>shapes = {"oval", "box", "polygon", "triangle", "egg"};

void label_leaf(
    const std::string& name,
    int idx,
    const std::vector<int>& cluster_assignments,
    const std::vector<int>& emissions,
    const std::string& indent,
    std::ofstream& s
) {
    const std::string& color = colors[cluster_assignments[idx] % colors.size()];
    const std::string& shape = shapes[emissions[idx] % shapes.size()];
    s << indent << name << " [label=\"" << idx << "\", fillcolor=" << color << ", shape=" << shape << "];\n";
}

void tree_to_dot(
    const char *file,
    const std::unordered_map<int, std::tuple<int, int>>& coal_tree,
    const std::vector<int>& cluster_assignments,
    const std::vector<int>& emissions,
    int l, int L
) {
    std::ofstream s = (l == 0) ? std::ofstream{file} : std::ofstream{file, std::ios::app};
    if (!s.is_open()) { throw std::runtime_error("Failed to open tree viz output file."); }

    if (l == 0) {
        s << "digraph G" << l << " {\n";
        s << "    node [style=filled]\n";
    }
    s << "    subgraph l" << l << " {" << '\n';

    const std::string indent(8, ' ');

    for (const auto& [node, children] : coal_tree) {
        const auto& [left, right] = children;

        std::string n_str = node_name(node, l);
        std::string l_str = node_name(left, l);
        std::string r_str  = node_name(right, l);

        s << indent << n_str << " [label=\"" << node << "\"];\n";
        s << indent << n_str << " -> " << l_str << ";\n";
        s << indent << n_str << " -> " << r_str << ";\n";

        if (left >= 0) { label_leaf(l_str, left, cluster_assignments, emissions, indent, s); }
        if (right >= 0) { label_leaf(r_str, right, cluster_assignments, emissions, indent, s); }
    }
    s << "    }\n";

    if (l == L-1) {
        s << "}\n";
    }
}

