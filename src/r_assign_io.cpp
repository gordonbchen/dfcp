#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <vector>
#include "r_assign_io.hpp"
#include "clusters.hpp"
#include "io.hpp"
#include "util.hpp"


namespace {

constexpr std::array<char, 4> r_assign_magic{'D', 'F', 'R', 'A'};
constexpr std::size_t r_assign_header_size = r_assign_magic.size() + 2 * sizeof(std::uint32_t);

}


RAssign read_r_assign(const char* fname) {
    std::ifstream input(fname, std::ios::binary);
    if (!input.is_open()) { throw std::runtime_error("Failed to open R assignment file."); }

    std::array<char, r_assign_magic.size()> magic{};
    read_bytes(input, magic.data(), magic.size());
    if (magic != r_assign_magic) { throw std::runtime_error("Invalid R assignment magic."); }

    std::array<std::uint32_t, 2> dimensions;
    read_bytes(input, dimensions.data(), sizeof(dimensions));
    auto [n_file, l_file] = dimensions;
    constexpr auto int_max = static_cast<std::uint32_t>(std::numeric_limits<int>::max());
    if (n_file == 0 || n_file > int_max || l_file == 0 || l_file > int_max) {
        throw std::runtime_error("R assignment dimensions must fit in positive int values.");
    }

    std::uintmax_t expected_size = r_assign_header_size
        + static_cast<std::uintmax_t>(n_file) * l_file * sizeof(std::uint32_t);
    if (std::filesystem::file_size(fname) != expected_size) {
        throw std::runtime_error("R assignment file size does not match its header.");
    }

    RAssign r_assign{
        .ids=std::vector<std::uint32_t>(static_cast<std::size_t>(n_file) * l_file),
        .N=static_cast<int>(n_file),
        .L=static_cast<int>(l_file),
    };
    read_bytes(input, r_assign.ids.data(), r_assign.ids.size() * sizeof(std::uint32_t));
    return r_assign;
}

void write_r_assign(const char* fname, const Clusters& clusters) {
    AtomicBinaryWriter output(fname);
    std::ostream& stream = output.stream();
    write_bytes(stream, r_assign_magic.data(), r_assign_magic.size());
    std::array<std::uint32_t, 2> dimensions{
        static_cast<std::uint32_t>(clusters.HP.N), static_cast<std::uint32_t>(clusters.HP.L)
    };
    write_bytes(stream, dimensions.data(), sizeof(dimensions));

    std::vector<std::uint32_t> row(clusters.HP.L);
    for (int i = 0; i < clusters.HP.N; ++i) {
        for (int l = 0; l < clusters.HP.L; ++l) {
            row[l] = clusters.r_assign[idx2d(i, l, clusters.HP.L)]->id;
        }
        write_bytes(stream, row.data(), row.size() * sizeof(std::uint32_t));
    }
    output.finish();
}
