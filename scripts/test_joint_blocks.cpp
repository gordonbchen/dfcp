#include <cassert>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "pbwt.hpp"
#include "seq_array.hpp"

int main() {
    // Every observed allele has a carrier, but target 1's two-locus pattern does not.
    SeqArray ref(3, 2);
    ref.x = {0b00, 0b01, 0b10};
    SeqArray target(2, 2);
    target.x = {0b00, 0b11};
    std::unordered_map<int, int> observed{{0, 0}, {1, 1}};

    auto individual = pbwt_blocks_greedy(ref, 4, &target, 0, &observed);
    auto joint = pbwt_blocks_greedy(ref, 4, &target, -1, &observed);
    assert((individual.snp_offsets == std::vector<int>{0, 2}));
    assert((joint.snp_offsets == std::vector<int>{0, 1, 2}));
    assert(joint.n_compat_cuts == 1);

    // A single observed allele absent from the reference cannot be accommodated.
    SeqArray absent_ref(2, 2);
    absent_ref.x = {0b00, 0b01};
    bool rejected = false;
    try {
        pbwt_blocks_greedy(absent_ref, 4, &target, -1, &observed);
    }
    catch (const std::runtime_error&) {
        rejected = true;
    }
    assert(rejected);
}
