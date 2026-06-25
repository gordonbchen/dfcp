#pragma once

#include <unordered_map>
#include <tuple>
#include <vector>


std::vector<int> parse_pos_file_idx(const char *fname, int start_idx, int L);

std::tuple<std::vector<std::unordered_map<int, std::tuple<int, int>>>, std::vector<int>> parse_tree_file(
    const char *fname
);

int calc_excess_parsimony(
    const std::unordered_map<int, std::tuple<int, int>>& coal_tree,
    const std::vector<int>& cluster_assignments,
    int n_clusters
);
