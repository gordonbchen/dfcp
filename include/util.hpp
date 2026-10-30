#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <utility>
#include "seq_array.hpp"


inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

inline int8_t get_xil(const SeqArray& x, int i, int l, const std::unordered_map<int, int> *obs_ls) {
    if (obs_ls == nullptr) {
        return x(i, l);
    }
    auto mapped_l = obs_ls->find(l);
    return mapped_l == obs_ls->end() ? -1 : x(i, mapped_l->second);
}


struct PairPointerHash {
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T*, U*>& p) const {
        auto hash1 = std::hash<T*>{}(p.first);
        auto hash2 = std::hash<U*>{}(p.second);
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
};


struct BitSet {
    size_t size;
    size_t n_words;
    std::vector<uint64_t> bm;
    size_t n_set;

    BitSet(size_t size_) : size(size_), n_words((size_ + 63) >> 6),
                           bm((size_ + 63) >> 6, 0), n_set(0) {}

    void set(size_t idx) {
        if (idx >= size) {
            throw std::invalid_argument("idx is too large for bit set.");
        }
        n_set += (bm[idx >> 6] & (1ULL << (((1 << 6) - 1) & idx))) == 0;
        bm[idx >> 6] |= 1ULL << (((1 << 6) - 1) & idx);
    }

    bool is_empty() {
        return n_set == 0;
    }

    void intersect(const BitSet& l, const BitSet& r) {
        if (!((size == l.size) && (l.size == r.size))) {
            throw std::invalid_argument("BitSets must be the same size.");
        }
        n_set = 0;
        for (size_t i = 0; i < n_words; ++i) {
            bm[i] = l.bm[i] & r.bm[i];
            n_set += std::popcount(bm[i]);
        }
    }

    void sunion(const BitSet& l, const BitSet& r) {
        if (!((size == l.size) && (l.size == r.size))) {
            throw std::invalid_argument("BitSets must be the same size.");
        }
        n_set = 0;
        for (size_t i = 0; i < n_words; ++i) {
            bm[i] = l.bm[i] | r.bm[i];
            n_set += std::popcount(bm[i]);
        }
    }

    bool operator==(const BitSet& other) const {
        if ((size != other.size) || (n_set != other.n_set)) { return false; }
        for (size_t i = 0; i < n_words; ++i) {
            if (bm[i] != other.bm[i]) {
                return false;
            }
        }
        return true;
    }
};

struct BitSetHash {
    size_t operator()(const BitSet& bs) const noexcept {
        size_t seed = std::hash<size_t>{}(bs.size);
        for (uint64_t word : bs.bm) {
            size_t h = std::hash<uint64_t>{}(word);
            seed ^= h + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


std::vector<size_t> count_emissions(const SeqArray& x, int K);
std::vector<int8_t> get_emission_modes(const std::vector<size_t>& emission_counts, const int L, const int K);
