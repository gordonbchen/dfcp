#pragma once

#include <vector>
#include <utility>
#include "seq_array.hpp"


std::pair<std::vector<int>, std::vector<int>> pbwt(const SeqArray& x);
