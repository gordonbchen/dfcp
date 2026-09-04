#pragma once

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iosfwd>


void read_bytes(std::istream& stream, void* data, std::size_t size);
void write_bytes(std::ostream& stream, const void* data, std::size_t size);


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
