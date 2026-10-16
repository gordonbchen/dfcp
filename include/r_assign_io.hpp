#pragma once

#include <cstdint>
#include <cstddef>
#include <filesystem>
#include <vector>
#include "clusters.hpp"


struct RAssign {
    std::vector<std::uint32_t> ids;
    int N;
    int L;

    std::uint32_t operator()(int i, int l) const {
        return ids[static_cast<size_t>(i) * L + l];
    }
};


RAssign read_r_assign(const std::filesystem::path& fname);
void write_r_assign(const std::filesystem::path& fname, const Clusters& clusters);
