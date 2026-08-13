#pragma once

#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>


std::vector<int> parse_pos_file_idx(const char *fname, int start_idx, int L);


struct CoalNode {
    int left;
    int right;
    double left_height;
    double right_height;
};

std::pair<std::vector<std::unordered_map<int, CoalNode>>, std::vector<int>> parse_tree_file(
    const char *fname
);


int count_observed_labels(const std::vector<int>& labels, size_t max_size);

int calc_excess_parsimony(
    const std::unordered_map<int, CoalNode>& coal_tree,
    const std::vector<int>& cluster_assignments,
    int n_clusters
);


std::vector<int> get_tree_idxs(const std::vector<int>& variant_pos, const std::vector<int>& recomb_pos);

std::pair<double, int> calc_max_clade_iou(
    const std::unordered_map<int, CoalNode>& coal_tree,
    const std::vector<int>& cluster_assign,
    int cluster_idx,
    int cluster_size
);

double calc_node_height(const std::unordered_map<int, CoalNode>& coal_tree, int idx);


void tree_to_dot(
    const char *file,
    const std::unordered_map<int, CoalNode>& coal_tree,
    const std::vector<int>& cluster_assignments,
    const std::vector<int>& emissions,
    int l,
    int L
);

