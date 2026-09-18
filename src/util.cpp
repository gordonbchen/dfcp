#include <cstdlib>
#include <stdexcept>
#include "util.hpp"


double parse_double(const char* value) {
    char* end = nullptr;
    double parsed = std::strtod(value, &end);
    if (end == value) { throw std::invalid_argument("Failed to parse double arg value."); }
    return parsed;
}

int parse_int(const char* value) {
    char* end = nullptr;
    int parsed = std::strtol(value, &end, 10);
    if (end == value) { throw std::invalid_argument("Failed to parse int arg value."); }
    return parsed;
}
