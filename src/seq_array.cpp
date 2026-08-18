#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <ios>
#include <vector>
#include "seq_array.hpp"


namespace {

constexpr std::uint32_t seq_file_magic = 0x50434644;  // "DFCP" as a little-endian integer.
constexpr int word_bits = 64;

std::uint32_t read_int32_le(std::ifstream& file) {
    std::uint8_t bytes[4];
    file.read(reinterpret_cast<char *>(bytes), sizeof(bytes));
    if (!file) {
        throw std::runtime_error("Truncated seq file header.");
    }
    return static_cast<std::uint32_t>(bytes[0])
        | (static_cast<std::uint32_t>(bytes[1]) << 8)
        | (static_cast<std::uint32_t>(bytes[2]) << 16)
        | (static_cast<std::uint32_t>(bytes[3]) << 24);
}

std::uint64_t byte_swap(std::uint64_t value) {
    value = ((value & 0x00ff00ff00ff00ffULL) << 8) | ((value >> 8) & 0x00ff00ff00ff00ffULL);
    value = ((value & 0x0000ffff0000ffffULL) << 16) | ((value >> 16) & 0x0000ffff0000ffffULL);
    return (value << 32) | (value >> 32);
}

void transpose_64(std::uint64_t *words) {
    // Transpose bits (row, column) to (column, row) using six block-swap stages.
    std::uint64_t mask = 0x00000000ffffffffULL;
    for (int shift = 32; shift != 0; shift >>= 1) {
        for (int i = 0; i < word_bits; i = (i + shift + 1) & ~shift) {
            std::uint64_t swap = ((words[i] >> shift) ^ words[i + shift]) & mask;
            words[i] ^= swap << shift;
            words[i + shift] ^= swap;
        }
        if (shift > 1) {
            mask ^= mask << (shift >> 1);
        }
    }
}

}


SeqArray::SeqArray(int N_, int L_) : N(N_), L(L_) {
    if (N <= 0 || L <= 0) {
        throw std::invalid_argument("SeqArray dimensions must be positive.");
    }
    x.resize(static_cast<std::size_t>(N) * words_for_loci(L), 0);
}

SeqArray read_seq_file(const char *filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open seq file.");
    }

    if (read_int32_le(file) != seq_file_magic) {
        throw std::runtime_error("Invalid seq file magic.");
    }

    std::uint32_t n_file = read_int32_le(file);
    std::uint32_t l_file = read_int32_le(file);
    constexpr auto int_max = static_cast<std::uint32_t>(std::numeric_limits<int>::max());
    if (n_file == 0 || n_file > int_max) {
        throw std::runtime_error("Invalid N in seq file.");
    }
    if (l_file == 0 || l_file > int_max) {
        throw std::runtime_error("Invalid L in seq file.");
    }

    int N = static_cast<int>(n_file);
    int L = static_cast<int>(l_file);
    SeqArray sequences(N, L);

    std::size_t words_per_locus = SeqArray::words_for_loci(N);
    std::size_t words_per_seq = sequences.words_per_seq();
    std::vector<std::uint64_t> locus_block(static_cast<std::size_t>(word_bits) * words_per_locus);

    for (int first_locus = 0; first_locus < L; first_locus += word_bits) {
        int block_loci = std::min(word_bits, L - first_locus);
        std::size_t block_words = static_cast<std::size_t>(block_loci) * words_per_locus;
        file.read(reinterpret_cast<char *>(locus_block.data()),
                  static_cast<std::streamsize>(block_words * sizeof(std::uint64_t)));
        if (!file) {
            throw std::runtime_error("Truncated bitpacked seq file.");
        }

        if constexpr (std::endian::native == std::endian::big) {
            for (std::size_t i = 0; i < block_words; ++i) {
                locus_block[i] = byte_swap(locus_block[i]);
            }
        }

        int used_bits = N % word_bits;
        if (used_bits != 0) {
            std::uint64_t padding_mask = ~((std::uint64_t{1} << used_bits) - 1);
            for (int l = 0; l < block_loci; ++l) {
                if ((locus_block[static_cast<std::size_t>(l) * words_per_locus
                                 + words_per_locus - 1] & padding_mask) != 0) {
                    throw std::runtime_error("Nonzero locus padding bits in seq file.");
                }
            }
        }

        for (std::size_t seq_word = 0; seq_word < words_per_locus; ++seq_word) {
            // File words are loci by sequences; each tile becomes sequences by loci.
            std::uint64_t tile[word_bits]{};
            for (int l = 0; l < block_loci; ++l) {
                tile[l] = locus_block[static_cast<std::size_t>(l) * words_per_locus + seq_word];
            }
            transpose_64(tile);

            int first_seq = static_cast<int>(seq_word * word_bits);
            int block_seqs = std::min(word_bits, N - first_seq);
            for (int i = 0; i < block_seqs; ++i) {
                sequences.x[static_cast<std::size_t>(first_seq + i) * words_per_seq
                            + first_locus / word_bits] = tile[i];
            }
        }
    }

    if (file.peek() != std::ifstream::traits_type::eof()) {
        throw std::runtime_error("Unexpected trailing data in seq file.");
    }
    return sequences;
}
