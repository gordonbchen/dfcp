#pragma once

#include <cstddef>
#include <vector>


size_t idx2d(size_t r, size_t c, size_t width);

std::vector<char> count_modes(const std::vector<char>& x, const int N, const int L, const int K);

