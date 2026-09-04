#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include "io.hpp"


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
