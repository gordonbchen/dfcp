#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "ref_filter.hpp"
#include "seq_array.hpp"


RefFilter::RefFilter(const SeqArray& ref_, const std::unordered_map<int, int>& obs_ls) :
    ref(ref_), observed_loci(obs_ls.size()),
    weight_zero(obs_ls.size()), weight_one(obs_ls.size())
{
    for (const auto& [locus, target_locus] : obs_ls) {
        observed_loci[target_locus] = locus;
    }

    for (int target_locus = 0; target_locus < static_cast<int>(observed_loci.size()); ++target_locus) {
        int n_one = 0;
        for (int i = 0; i < ref.N; ++i) {
            n_one += ref(i, observed_loci[target_locus]);
        }
        int n_zero = ref.N - n_one;
        // An absent allele cannot match any reference, so it never needs a log weight.
        if (n_zero != 0) {
            weight_zero[target_locus] = -std::log(static_cast<double>(n_zero) / ref.N);
        }
        if (n_one != 0) {
            weight_one[target_locus] = -std::log(static_cast<double>(n_one) / ref.N);
        }
    }
}


std::vector<int> RefFilter::top_references(
    const SeqArray& target, int target_idx, int top_k
) const {
    if (top_k <= 0 || target_idx < 0 || target_idx >= target.N
        || target.L != static_cast<int>(observed_loci.size())) {
        throw std::invalid_argument("Invalid target or reference top K for filtering.");
    }

    std::vector<int> indexes(ref.N);
    std::iota(indexes.begin(), indexes.end(), 0);
    if (top_k >= ref.N) {
        return indexes;
    }

    std::vector<double> scores(ref.N, 0.0);
    for (int i = 0; i < ref.N; ++i) {
        double score = 0.0;
        for (int target_locus = 0; target_locus < target.L; ++target_locus) {
            int allele = target(target_idx, target_locus);
            if (ref(i, observed_loci[target_locus]) == allele) {
                score += allele == 0 ? weight_zero[target_locus] : weight_one[target_locus];
            }
        }
        scores[i] = score;
    }

    std::partial_sort(indexes.begin(), indexes.begin() + top_k, indexes.end(),
                      [&](int a, int b) {
                          return scores[a] > scores[b] || (scores[a] == scores[b] && a < b);
                      });
    indexes.resize(top_k);

    // Every observed target allele needs at least one carrier for block construction.
    std::vector<bool> selected(ref.N, false);
    for (int i : indexes) {
        selected[i] = true;
    }
    for (int target_locus = 0; target_locus < target.L; ++target_locus) {
        int allele = target(target_idx, target_locus);
        int locus = observed_loci[target_locus];
        bool covered = false;
        for (int i : indexes) {
            if (ref(i, locus) == allele) {
                covered = true;
                break;
            }
        }
        if (covered) {
            continue;
        }

        int best = -1;
        for (int i = 0; i < ref.N; ++i) {
            if (!selected[i] && ref(i, locus) == allele
                && (best == -1 || scores[i] > scores[best]
                    || (scores[i] == scores[best] && i < best))) {
                best = i;
            }
        }
        if (best == -1) {
            throw std::runtime_error("An observed target allele is absent from the full reference panel.");
        }
        indexes.push_back(best);
        selected[best] = true;
    }

    std::sort(indexes.begin(), indexes.end());
    return indexes;
}


SeqArray select_references(const SeqArray& ref, const std::vector<int>& indexes) {
    if (indexes.empty()) {
        throw std::invalid_argument("The selected reference panel is empty.");
    }
    SeqArray selected(static_cast<int>(indexes.size()), ref.L);
    std::size_t words_per_seq = ref.words_per_seq();
    for (int i = 0; i < selected.N; ++i) {
        int source = indexes[i];
        if (source < 0 || source >= ref.N) {
            throw std::invalid_argument("Selected reference index is out of range.");
        }
        std::copy_n(ref.x.begin() + static_cast<std::size_t>(source) * words_per_seq,
                    words_per_seq,
                    selected.x.begin() + static_cast<std::size_t>(i) * words_per_seq);
    }
    return selected;
}
