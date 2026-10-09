#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "ref_filter.hpp"
#include "seq_array.hpp"


int main() {
    SeqArray ref(5, 5);
    ref.x = {0b0000, 0b0001, 0b0010, 0b0011, 0b1111};
    SeqArray target(2, 5);
    target.x = {0b0001, 0b0100};
    // Locus 4 is monomorphic zero, so its matches carry zero weight.
    std::unordered_map<int, int> observed{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};

    RefFilter filter(ref, observed);
    assert((filter.top_references(target, 0, 2) == std::vector<int>{0, 1}));
    // The top-scoring reference alone lacks target 1 at locus 2; add its carrier.
    assert((filter.top_references(target, 1, 1) == std::vector<int>{0, 4}));
    assert((filter.top_references(target, 0, 5) == std::vector<int>{0, 1, 2, 3, 4}));

    SeqArray selected = select_references(ref, {0, 4});
    assert(selected.N == 2 && selected.L == 5);
    assert(selected.x[0] == ref.x[0] && selected.x[1] == ref.x[4]);
}
