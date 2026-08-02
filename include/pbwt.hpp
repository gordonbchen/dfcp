#pragma once

#include <vector>
#include <utility>
#include <cstdint>


std::pair<std::vector<int>, std::vector<int>> pbwt(const std::vector<int8_t>& x, int N, int L);
std::pair<std::vector<int>, std::vector<int>> reverse_pbwt(const std::vector<int8_t>& x, int N, int L);

