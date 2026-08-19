#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <istream>
#include <ostream>
#include <span>
#include <stdexcept>
#include <system_error>
#include <vector>
#include "io.hpp"


namespace {

std::uint64_t byte_swap(std::uint64_t value) {
    value = ((value & 0x00ff00ff00ff00ffULL) << 8) | ((value >> 8) & 0x00ff00ff00ff00ffULL);
    value = ((value & 0x0000ffff0000ffffULL) << 16) | ((value >> 16) & 0x0000ffff0000ffffULL);
    return (value << 32) | (value >> 32);
}

}


void read_bytes(std::istream& stream, void* data, std::size_t size) {
    stream.read(static_cast<char*>(data), static_cast<std::streamsize>(size));
    if (!stream) {
        throw std::runtime_error("Unexpected end of binary file.");
    }
}

void write_bytes(std::ostream& stream, const void* data, std::size_t size) {
    stream.write(static_cast<const char*>(data), static_cast<std::streamsize>(size));
    if (!stream) {
        throw std::runtime_error("Failed to write binary file.");
    }
}


std::uint32_t read_uint32_le(std::istream& stream) {
    std::array<std::uint8_t, 4> bytes{};
    read_bytes(stream, bytes.data(), bytes.size());
    return static_cast<std::uint32_t>(bytes[0])
        | (static_cast<std::uint32_t>(bytes[1]) << 8)
        | (static_cast<std::uint32_t>(bytes[2]) << 16)
        | (static_cast<std::uint32_t>(bytes[3]) << 24);
}

void write_uint32_le(std::ostream& stream, std::uint32_t value) {
    std::array<std::uint8_t, 4> bytes{
        static_cast<std::uint8_t>(value),
        static_cast<std::uint8_t>(value >> 8),
        static_cast<std::uint8_t>(value >> 16),
        static_cast<std::uint8_t>(value >> 24),
    };
    write_bytes(stream, bytes.data(), bytes.size());
}


void read_uint16s_le(std::istream& stream, std::span<std::uint16_t> values) {
    read_bytes(stream, values.data(), values.size_bytes());
    if constexpr (std::endian::native == std::endian::big) {
        for (std::uint16_t& value : values) {
            value = static_cast<std::uint16_t>((value << 8) | (value >> 8));
        }
    }
}

void write_uint16s_le(std::ostream& stream, std::span<const std::uint16_t> values) {
    if constexpr (std::endian::native == std::endian::little) {
        write_bytes(stream, values.data(), values.size_bytes());
    }
    else {
        std::vector<std::uint16_t> swapped(values.begin(), values.end());
        for (std::uint16_t& value : swapped) {
            value = static_cast<std::uint16_t>((value << 8) | (value >> 8));
        }
        write_bytes(stream, swapped.data(), swapped.size() * sizeof(std::uint16_t));
    }
}


void read_uint64s_le(std::istream& stream, std::span<std::uint64_t> values) {
    read_bytes(stream, values.data(), values.size_bytes());
    if constexpr (std::endian::native == std::endian::big) {
        for (std::uint64_t& value : values) {
            value = byte_swap(value);
        }
    }
}

void write_float32s_le(std::ostream& stream, std::span<const float> values) {
    std::vector<std::uint32_t> bits(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        bits[i] = std::bit_cast<std::uint32_t>(values[i]);
    }
    if constexpr (std::endian::native == std::endian::big) {
        for (std::uint32_t& value : bits) {
            value = (value << 24) | ((value << 8) & 0x00ff0000) | ((value >> 8) & 0x0000ff00)
                | (value >> 24);
        }
    }
    write_bytes(stream, bits.data(), bits.size() * sizeof(std::uint32_t));
}


AtomicBinaryWriter::AtomicBinaryWriter(const char* filename) :
    output_path(filename), temporary_path(output_path.string() + ".tmp"),
    output(temporary_path, std::ios::binary | std::ios::trunc)
{
    if (!output.is_open()) {
        throw std::runtime_error("Failed to open temporary output file.");
    }
}

AtomicBinaryWriter::~AtomicBinaryWriter() {
    if (!finished) {
        output.close();
        std::error_code error;
        std::filesystem::remove(temporary_path, error);
    }
}

std::ostream& AtomicBinaryWriter::stream() {
    return output;
}

void AtomicBinaryWriter::finish() {
    output.close();
    if (!output) {
        throw std::runtime_error("Failed to close output file.");
    }
    std::error_code error;
    std::filesystem::rename(temporary_path, output_path, error);
    if (error) {
        throw std::runtime_error("Failed to rename output file: " + error.message());
    }
    finished = true;
}
