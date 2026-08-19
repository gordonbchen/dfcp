#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iosfwd>
#include <span>


void read_bytes(std::istream& stream, void* data, std::size_t size);
void write_bytes(std::ostream& stream, const void* data, std::size_t size);

std::uint32_t read_uint32_le(std::istream& stream);
void write_uint32_le(std::ostream& stream, std::uint32_t value);
void read_uint16s_le(std::istream& stream, std::span<std::uint16_t> values);
void write_uint16s_le(std::ostream& stream, std::span<const std::uint16_t> values);
void read_uint64s_le(std::istream& stream, std::span<std::uint64_t> values);
void write_float32s_le(std::ostream& stream, std::span<const float> values);


class AtomicBinaryWriter {
    std::filesystem::path output_path;
    std::filesystem::path temporary_path;
    std::ofstream output;
    bool finished = false;

    public:
        explicit AtomicBinaryWriter(const char* filename);
        ~AtomicBinaryWriter();

        AtomicBinaryWriter(const AtomicBinaryWriter&) = delete;
        AtomicBinaryWriter& operator=(const AtomicBinaryWriter&) = delete;

        std::ostream& stream();
        void finish();
};
