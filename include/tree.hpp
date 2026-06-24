#pragma once

#include <unordered_map>
#include <tuple>
#include <vector>
#include "clusters.hpp"


std::vector<std::unordered_map<int, std::tuple<int, int>>> parse_tree_file(const char *tree_file_name, int L);

int calc_excess_parsimony(
    int l,
    const std::unordered_map<int, std::tuple<int, int>>& coal_tree,
    const Clusters& clusters,
    const std::unordered_map<Cluster*, int>& cluster_idxs
);

