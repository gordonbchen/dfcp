#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>


struct SeqArray {
    std::vector<std::uint64_t> x;
    int N;
    int L;

    SeqArray(int N_, int L_);

    static std::size_t words_for_loci(int n_loci) {
        return (static_cast<std::size_t>(n_loci) + 63) / 64;
    }

    std::size_t words_per_seq() const {
        return words_for_loci(L);
    }

    std::int8_t operator()(int i, int l) const {
        std::size_t word_idx = static_cast<std::size_t>(i) * words_per_seq() + l / 64;
        return static_cast<std::int8_t>((x[word_idx] >> (l % 64)) & 1);
    }
};


SeqArray read_seq_file(const char *filename);
