#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>
#include "io.hpp"
#include "seq_array.hpp"


namespace {

constexpr std::array<char, 4> seq_file_magic{'D', 'F', 'C', 'P'};
constexpr std::size_t seq_header_size = seq_file_magic.size() + 2 * sizeof(std::uint32_t);
constexpr int word_bits = 64;

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


SeqArray::SeqArray(int N_, int L_) :
    x(static_cast<std::size_t>(N_) * words_for_loci(L_), 0), N(N_), L(L_)
{}

SeqArray read_seq_file(const char *filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open seq file.");
    }

    std::array<char, seq_file_magic.size()> magic{};
    read_bytes(file, magic.data(), magic.size());
    if (magic != seq_file_magic) {
        throw std::runtime_error("Invalid seq file magic.");
    }

    std::uint32_t n_file = read_uint32_le(file);
    std::uint32_t l_file = read_uint32_le(file);
    constexpr auto int_max = static_cast<std::uint32_t>(std::numeric_limits<int>::max());
    if (n_file == 0 || n_file > int_max || l_file == 0 || l_file > int_max) {
        throw std::runtime_error("Seq file dimensions must fit in positive int values.");
    }

    int N = static_cast<int>(n_file);
    int L = static_cast<int>(l_file);
    SeqArray sequences(N, L);

    std::size_t words_per_locus = SeqArray::words_for_loci(N);
    std::size_t words_per_seq = sequences.words_per_seq();
    std::uint64_t n_file_words = static_cast<std::uint64_t>(L) * words_per_locus;
    std::uintmax_t file_size = std::filesystem::file_size(filename);
    if (file_size < seq_header_size || (file_size - seq_header_size) % sizeof(std::uint64_t) != 0
        || (file_size - seq_header_size) / sizeof(std::uint64_t) != n_file_words) {
        throw std::runtime_error("Seq file size does not match its header.");
    }

    std::vector<std::uint64_t> locus_block(static_cast<std::size_t>(word_bits) * words_per_locus);

    for (int first_locus = 0; first_locus < L; first_locus += word_bits) {
        int block_loci = std::min(word_bits, L - first_locus);
        std::size_t block_words = static_cast<std::size_t>(block_loci) * words_per_locus;
        read_uint64s_le(file, std::span(locus_block.data(), block_words));

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
    return sequences;
}
