#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ios>
#include <stdexcept>
#include <string>
#include <unordered_map>
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

const char* parse_coal_subtree(const char *s, std::unordered_map<int, CoalNode>& coal_tree, int idx) {
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

    while (*s != ':') { ++s; }
    ++s;
    double left_height = std::strtod(s, nullptr);

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

    while (*s != ':') { ++s; }
    ++s;
    double right_height = std::strtod(s, nullptr);

    coal_tree.emplace(idx, CoalNode{left_idx, right_idx, left_height, right_height});
    return s;
}

int parse_coal_tree(const char *s, std::unordered_map<int, CoalNode>& coal_tree) {
    s = std::strstr(s, "pos_");
    if (s == nullptr) { throw std::runtime_error("Failed to parse tree, could not find 'pos_'"); }
    s += 4;

    int recomb_pos = std::strtol(s, nullptr, 10);
    parse_coal_subtree(s, coal_tree, -1);
    return recomb_pos;
}

std::pair<std::vector<std::unordered_map<int, CoalNode>>, std::vector<int>> parse_tree_file(
    const char *fname
) {
    std::ifstream tree_file{fname};
    if (!tree_file.is_open()) { throw std::runtime_error("Failed to open tree file."); };

    std::string line;
    for (int i = 0; i < 3; ++i) {
        std::getline(tree_file, line);
    }

    std::vector<std::unordered_map<int, CoalNode>> trees;
    std::vector<int> recomb_pos;
    int l = 0;
    while (std::getline(tree_file, line)) {
        if (line[0] == 'e') { break; }
        trees.emplace_back();
        recomb_pos.push_back(parse_coal_tree(line.c_str(), trees[l]));
        ++l;
    }
    return {std::move(trees), std::move(recomb_pos)};
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
    const std::unordered_map<int, CoalNode>& coal_tree,
    const std::vector<int>& cluster_assignments,
    const size_t n_clusters
) {
    BitSet cluster_bm{n_clusters};
    if (!coal_tree.contains(idx)) {
        cluster_bm.set(cluster_assignments[idx]);
        return ParsimonyMsg{0, std::move(cluster_bm)};
    }

    const auto& coal_node = coal_tree.at(idx);
    ParsimonyMsg lp = calc_parsimony(coal_node.left, coal_tree, cluster_assignments, n_clusters);
    ParsimonyMsg rp = calc_parsimony(coal_node.right, coal_tree, cluster_assignments, n_clusters);

    cluster_bm.intersect(lp.cluster_bm, rp.cluster_bm);
    if (!cluster_bm.is_empty()) {
        return ParsimonyMsg{lp.score + rp.score, std::move(cluster_bm)};
    }
    cluster_bm.sunion(lp.cluster_bm, rp.cluster_bm);
    return ParsimonyMsg{lp.score + rp.score + 1, std::move(cluster_bm)};
}

int calc_excess_parsimony(
    const std::unordered_map<int, CoalNode>& coal_tree,
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


struct CladeIOUMsg {
    double iou;
    int leaves;
    int isect;
    int root;
};

CladeIOUMsg calc_max_clade_iou_dfs(
    int v,
    const std::unordered_map<int, CoalNode>& coal_tree,
    const std::vector<int>& cluster_assign,
    int cluster_idx,
    int cluster_size
) {
    if (v < 0) {
        const auto& coal_node = coal_tree.at(v);
        CladeIOUMsg l_msg{calc_max_clade_iou_dfs(coal_node.left, coal_tree, cluster_assign, cluster_idx, cluster_size)};
        CladeIOUMsg r_msg{calc_max_clade_iou_dfs(coal_node.right, coal_tree, cluster_assign, cluster_idx, cluster_size)};

        CladeIOUMsg msg;
        msg.leaves = l_msg.leaves + r_msg.leaves;
        msg.isect = l_msg.isect + r_msg.isect;

        double iou = static_cast<double>(msg.isect) / (msg.leaves + cluster_size - msg.isect);
        msg.iou = std::max(iou, std::max(l_msg.iou, r_msg.iou));
        msg.root = (iou > std::max(l_msg.iou, r_msg.iou)) ? v : ((l_msg.iou > r_msg.iou) ? l_msg.root : r_msg.root);
        return msg;
    }

    CladeIOUMsg msg;
    msg.leaves = 1;
    msg.isect = cluster_assign[v] == cluster_idx ? 1 : 0;
    msg.iou = static_cast<double>(msg.isect) / (msg.leaves + cluster_size - msg.isect);
    msg.root = v;
    return msg;
}

std::pair<double, int> calc_max_clade_iou(
    const std::unordered_map<int, CoalNode>& coal_tree,
    const std::vector<int>& cluster_assign,
    int cluster_idx,
    int cluster_size
) {
    CladeIOUMsg msg = calc_max_clade_iou_dfs(-1, coal_tree, cluster_assign, cluster_idx, cluster_size);
    return {msg.iou, msg.root};
}

double calc_node_height(const std::unordered_map<int, CoalNode>& coal_tree, int idx) {
    if (idx >= 0) {
        return 0.0;
    }
    const CoalNode& node = coal_tree.at(idx);
    return node.left_height + calc_node_height(coal_tree, node.left);
}


std::string node_name(int x, int l) {
    return std::format("\"l{}_{}\"", l, x);
}

struct RGB {
    double r;
    double g;
    double b;
};

double hue_to_rgb(double p, double q, double t) {
    if (t < 0.0) {
        t += 1.0;
    }
    if (t > 1.0) {
        t -= 1.0;
    }

    if (t < 1.0 / 6.0) {
        return p + (q - p) * 6.0 * t;
    }
    if (t < 1.0 / 2.0) {
        return q;
    }
    if (t < 2.0 / 3.0) {
        return p + (q - p) * (2.0 / 3.0 - t) * 6.0;
    }

    return p;
}

RGB hsl_to_rgb(double h, double s, double l) {
    if (s == 0.0) {
        return RGB{l, l, l};
    }

    double q = (l < 0.5) ? l * (1.0 + s) : l + s - l * s;
    double p = 2.0 * l - q;

    double r = hue_to_rgb(p, q, h + 1.0 / 3.0);
    double g = hue_to_rgb(p, q, h);
    double b = hue_to_rgb(p, q, h - 1.0 / 3.0);

    return RGB{r, g, b};
}

std::string rgb_to_hex(RGB c) {
    int r = static_cast<int>(std::round(255.0 * std::clamp(c.r, 0.0, 1.0)));
    int g = static_cast<int>(std::round(255.0 * std::clamp(c.g, 0.0, 1.0)));
    int b = static_cast<int>(std::round(255.0 * std::clamp(c.b, 0.0, 1.0)));

    return std::format("#{:02X}{:02X}{:02X}", r, g, b);
}

std::string cluster_color(size_t cluster_idx) {
    constexpr double golden_ratio_conjugate = 0.6180339887498948482;
    double h = std::fmod(0.17 + golden_ratio_conjugate * cluster_idx, 1.0);
    double s = 0.7;
    double l = 0.7;

    return rgb_to_hex(hsl_to_rgb(h, s, l));
}

const std::vector<std::string>shapes = {"oval", "box", "polygon", "triangle", "egg"};

void label_leaf(
    const std::string& name,
    int idx,
    const std::vector<int>& cluster_assignments,
    const std::vector<int>& emissions,
    const std::string& indent,
    std::ofstream& s
) {
    const std::string& color = cluster_color(cluster_assignments[idx]);
    const std::string& shape = shapes[emissions[idx] % shapes.size()];
    s << indent << name << " [label=\"" << idx << "\", fillcolor=\"" << color << "\", shape=" << shape << "];\n";
}

void tree_to_dot(
    const char *file,
    const std::unordered_map<int, CoalNode>& coal_tree,
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

    for (const auto& [node_idx, coal_node] : coal_tree) {
        std::string n_str = node_name(node_idx, l);
        std::string l_str = node_name(coal_node.left, l);
        std::string r_str  = node_name(coal_node.right, l);

        s << indent << n_str << " [label=\"" << node_idx << "\"];\n";
        s << indent << n_str << " -> " << l_str << ";\n";
        s << indent << n_str << " -> " << r_str << ";\n";

        if (coal_node.left >= 0) { label_leaf(l_str, coal_node.left, cluster_assignments, emissions, indent, s); }
        if (coal_node.right >= 0) { label_leaf(r_str, coal_node.right, cluster_assignments, emissions, indent, s); }
    }
    s << "    }\n";

    if (l == L-1) {
        s << "}\n";
    }
}
