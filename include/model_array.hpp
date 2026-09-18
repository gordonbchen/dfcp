#pragma once

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <unordered_map>

#include "pbwt.hpp"
#include "seq_array.hpp"

struct ModelArray {
    const SeqArray& snps;
    std::optional<PbwtBlocks> blocks;
    int L;

    ModelArray(
        const SeqArray& snps_, int block_max_k_,
        const SeqArray* target,
        const std::unordered_map<int, int>* obs_ls
    ) :
        snps(snps_), L(snps_.L)
    {
        if (block_max_k_ < 0) {
            throw std::invalid_argument("block_max_k must be nonnegative.");
        }
        if (block_max_k_ > 0) {
            blocks = pbwt_blocks_greedy(snps, block_max_k_, target, obs_ls);
            L = blocks->L;
        }
    }

    std::uint32_t operator()(int i, int l) const {
        return blocks.has_value() ? (*blocks)(i, l) : snps(i, l);
    }

    bool is_blocked() const {
        return blocks.has_value();
    }

    int n_emissions(int l) const {
        return blocks.has_value() ? blocks->n_emissions(l) : 2;
    }

    int snp_start(int l) const {
        return blocks.has_value() ? blocks->snp_start(l) : l;
    }

    int snp_end(int l) const {
        return blocks.has_value() ? blocks->snp_end(l) : l + 1;
    }

    int rep(int l, int emission) const {
        return blocks->rep(l, emission);
    }
};
