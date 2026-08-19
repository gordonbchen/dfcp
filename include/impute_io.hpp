#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <span>
#include <vector>
#include "io.hpp"


struct ImputeProbHeader {
    std::uint32_t n_sequences;
    std::uint32_t n_loci;
};


class ImputeProbReader {
    static constexpr std::size_t header_size = 4 + 2 * sizeof(std::uint32_t);

    std::ifstream input;
    ImputeProbHeader file_header;

    public:
        explicit ImputeProbReader(const char* filename);

        const ImputeProbHeader& header() const;
        void read_row(std::span<std::uint16_t> values);
};


class ImputeProbWriter {
    AtomicBinaryWriter output;
    std::uint32_t n_sequences;
    std::uint32_t n_loci;
    std::uint32_t rows_written = 0;
    std::vector<std::uint16_t> row_buffer;

    public:
        ImputeProbWriter(const char* filename, int n_sequences_, int n_loci_);

        ImputeProbWriter(const ImputeProbWriter&) = delete;
        ImputeProbWriter& operator=(const ImputeProbWriter&) = delete;

        void write_row(std::span<const double> allele_probs);
        void finish();
};


void write_impute_eval_file(
    const char* filename, std::span<const float> r2, std::span<const float> accuracy
);
