#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <span>
#include <stdexcept>
#include "impute_io.hpp"
#include "io.hpp"


namespace {

constexpr std::array<char, 4> prob_magic{'D', 'F', 'I', 'P'};
constexpr std::uint32_t fixed_point_max = std::numeric_limits<std::uint16_t>::max();

}


ImputeProbReader::ImputeProbReader(const char* filename) : input(filename, std::ios::binary) {
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open imputation probability file.");
    }

    std::array<char, prob_magic.size()> magic{};
    read_bytes(input, magic.data(), magic.size());
    if (magic != prob_magic) {
        throw std::runtime_error("Invalid imputation probability magic.");
    }

    std::array<std::uint32_t, 2> dimensions;
    read_bytes(input, dimensions.data(), sizeof(dimensions));
    file_header = {dimensions[0], dimensions[1]};
    if (file_header.n_sequences == 0 || file_header.n_loci == 0) {
        throw std::runtime_error("Imputation probability dimensions must be positive.");
    }

    std::uintmax_t file_size = std::filesystem::file_size(filename);
    std::uint64_t n_values = static_cast<std::uint64_t>(file_header.n_sequences) * file_header.n_loci;
    if (file_size < header_size || (file_size - header_size) % sizeof(std::uint16_t) != 0
        || (file_size - header_size) / sizeof(std::uint16_t) != n_values) {
        throw std::runtime_error("Imputation probability file size does not match its header.");
    }
}

const ImputeProbHeader& ImputeProbReader::header() const {
    return file_header;
}

void ImputeProbReader::read_row(std::span<std::uint16_t> values) {
    read_bytes(input, values.data(), values.size_bytes());
}


ImputeProbWriter::ImputeProbWriter(const char* filename, int n_sequences_, int n_loci_) :
    output(filename), n_sequences(n_sequences_), n_loci(n_loci_), row_buffer(n_loci)
{
    write_bytes(output.stream(), prob_magic.data(), prob_magic.size());
    std::array<std::uint32_t, 2> dimensions{n_sequences, n_loci};
    write_bytes(output.stream(), dimensions.data(), sizeof(dimensions));
}

void ImputeProbWriter::write_row(std::span<const double> allele_probs) {
    for (std::size_t l = 0; l < n_loci; ++l) {
        double probability = allele_probs[2 * l + 1];
        if (!std::isfinite(probability) || probability < 0.0 || probability > 1.0) {
            throw std::runtime_error("Imputation probability is outside [0, 1].");
        }
        row_buffer[l] = static_cast<std::uint16_t>(std::lround(probability * fixed_point_max));
    }
    write_bytes(output.stream(), row_buffer.data(), row_buffer.size() * sizeof(std::uint16_t));
    ++rows_written;
}

void ImputeProbWriter::finish() {
    if (rows_written != n_sequences) {
        throw std::runtime_error("Imputation probability file has the wrong number of rows.");
    }
    output.finish();
}
