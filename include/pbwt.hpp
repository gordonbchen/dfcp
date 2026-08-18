#pragma once

#include <vector>
#include <utility>
#include "seq_array.hpp"


std::pair<std::vector<int>, std::vector<int>> pbwt(const SeqArray& x);
std::pair<std::vector<int>, std::vector<int>> reverse_pbwt(const SeqArray& x);
