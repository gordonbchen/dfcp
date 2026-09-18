#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "pbwt.hpp"
#include "seq_array.hpp"
#include "util.hpp"

namespace {

struct PbwtStep {
    int n_patterns;
    int n_alleles;
};

struct StreamingPbwt {
    const SeqArray& x;
    std::vector<int> a;
    std::vector<int> d;
    std::vector<int> next_a;
    std::vector<int> next_d;
    std::vector<int> zeros;
    std::vector<int> ones;
    std::vector<int> zeros_d;
    std::vector<int> ones_d;

    explicit StreamingPbwt(const SeqArray& x_) :
        x(x_), a(x.N), d(x.N, 0), next_a(x.N), next_d(x.N),
        zeros(x.N), ones(x.N), zeros_d(x.N), ones_d(x.N)
    {
        for (int i = 0; i < x.N; ++i) {
            a[i] = i;
        }
    }

    PbwtStep advance(int snp, int block_start) {
        int n_zeros = 0;
        int n_ones = 0;
        int p = snp + 1;
        int q = snp + 1;
        for (int j = 0; j < x.N; ++j) {
            p = std::max(p, d[j]);
            q = std::max(q, d[j]);
            int seq = a[j];
            if (x(seq, snp) == 0) {
                zeros[n_zeros] = seq;
                zeros_d[n_zeros] = p;
                p = 0;
                ++n_zeros;
            }
            else {
                ones[n_ones] = seq;
                ones_d[n_ones] = q;
                q = 0;
                ++n_ones;
            }
        }

        int n_patterns = 1;
        for (int j = 0; j < n_zeros; ++j) {
            next_a[j] = zeros[j];
            next_d[j] = zeros_d[j];
            if (j > 0 && next_d[j] > block_start) {
                ++n_patterns;
            }
        }
        for (int j = 0; j < n_ones; ++j) {
            int next_j = n_zeros + j;
            next_a[next_j] = ones[j];
            next_d[next_j] = ones_d[j];
            if (next_j > 0 && next_d[next_j] > block_start) {
                ++n_patterns;
            }
        }
        a.swap(next_a);
        d.swap(next_d);
        return PbwtStep{n_patterns, (n_zeros != 0) + (n_ones != 0)};
    }

    void advance(int snp) {
        advance(snp, 0);
    }
};

class TargetCompat {
    const SeqArray* target;
    const std::unordered_map<int, int>* obs_ls;
    int n_ref;
    int n_ref_words;
    std::vector<std::uint64_t> candidates;
    std::vector<std::uint64_t> ref_ones;

    public:
        TargetCompat(
            int n_ref_, const SeqArray* target_,
            const std::unordered_map<int, int>* obs_ls_
        ) :
            target(target_),
            obs_ls(obs_ls_),
            n_ref(n_ref_),
            n_ref_words((n_ref_ + 63) / 64),
            candidates(target == nullptr ? 0 : static_cast<std::size_t>(target->N) * n_ref_words),
            ref_ones(target == nullptr ? 0 : n_ref_words)
        {
            if ((target == nullptr) != (obs_ls == nullptr)) {
                throw std::invalid_argument("target and obs_ls must either be both present or absent.");
            }
            reset();
        }

        void reset() {
            if (target == nullptr) {
                return;
            }
            std::fill(candidates.begin(), candidates.end(), ~std::uint64_t{0});
            int tail_bits = n_ref % 64;
            if (tail_bits == 0) {
                return;
            }
            std::uint64_t tail_mask = (std::uint64_t{1} << tail_bits) - 1;
            for (int i = 0; i < target->N; ++i) {
                candidates[static_cast<std::size_t>(i + 1) * n_ref_words - 1] = tail_mask;
            }
        }

        bool constrain(const SeqArray& ref, int snp) {
            if (target == nullptr) {
                return true;
            }
            auto obs_it = obs_ls->find(snp);
            if (obs_it == obs_ls->end()) {
                return true;
            }

            std::fill(ref_ones.begin(), ref_ones.end(), 0);
            for (int i = 0; i < ref.N; ++i) {
                if (ref(i, snp) != 0) {
                    ref_ones[i / 64] |= std::uint64_t{1} << (i % 64);
                }
            }

            for (int i = 0; i < target->N; ++i) {
                bool any = false;
                bool allele = (*target)(i, obs_it->second) != 0;
                std::size_t row = static_cast<std::size_t>(i) * n_ref_words;
                for (int w = 0; w < n_ref_words; ++w) {
                    std::uint64_t allele_mask = allele ? ref_ones[w] : ~ref_ones[w];
                    candidates[row + w] &= allele_mask;
                    any = any || candidates[row + w] != 0;
                }
                if (!any) {
                    return false;
                }
            }
            return true;
        }
};

PbwtBlocks label_blocks(
    const SeqArray& x, std::vector<int> snp_offsets,
    int n_max_k_cuts, int n_compat_cuts
) {
    int n_blocks = static_cast<int>(snp_offsets.size()) - 1;
    PbwtBlocks blocks{
        .labels=std::vector<std::uint32_t>(static_cast<std::size_t>(x.N) * n_blocks),
        .rep_offsets={0},
        .reps={},
        .emission_sizes={},
        .snp_offsets=std::move(snp_offsets),
        .L=n_blocks,
        .n_max_k_cuts=n_max_k_cuts,
        .n_compat_cuts=n_compat_cuts,
    };

    StreamingPbwt state{x};
    int block = 0;
    for (int snp = 0; snp < x.L; ++snp) {
        state.advance(snp);
        if (snp + 1 != blocks.snp_end(block)) {
            continue;
        }

        std::uint32_t emission = 0;
        int group_start = 0;
        blocks.reps.push_back(state.a[0]);
        blocks.labels[idx2d(state.a[0], block, n_blocks)] = emission;
        for (int j = 1; j < x.N; ++j) {
            if (state.d[j] > blocks.snp_start(block)) {
                blocks.emission_sizes.push_back(j - group_start);
                group_start = j;
                ++emission;
                blocks.reps.push_back(state.a[j]);
            }
            blocks.labels[idx2d(state.a[j], block, n_blocks)] = emission;
        }
        blocks.emission_sizes.push_back(x.N - group_start);
        blocks.rep_offsets.push_back(blocks.reps.size());
        ++block;
    }
    return blocks;
}

}

std::pair<std::vector<int>, std::vector<int>> pbwt(const SeqArray& x) {
    std::vector<int> a(static_cast<std::size_t>(x.N) * (x.L + 1));
    std::vector<int> d(static_cast<std::size_t>(x.N) * (x.L + 1), 0);
    StreamingPbwt state{x};
    for (int i = 0; i < x.N; ++i) {
        a[idx2d(i, 0, x.L + 1)] = state.a[i];
    }
    for (int l = 0; l < x.L; ++l) {
        state.advance(l);
        for (int i = 0; i < x.N; ++i) {
            a[idx2d(i, l + 1, x.L + 1)] = state.a[i];
            d[idx2d(i, l + 1, x.L + 1)] = state.d[i];
        }
    }
    return {std::move(a), std::move(d)};
}

PbwtBlocks pbwt_blocks_greedy(
    const SeqArray& ref, int max_k,
    const SeqArray* target,
    const std::unordered_map<int, int>* obs_ls
) {
    if (max_k <= 0) {
        throw std::invalid_argument("block_max_k must be positive.");
    }

    StreamingPbwt state{ref};
    TargetCompat target_compat{ref.N, target, obs_ls};
    std::vector<int> snp_offsets{0};
    int block_start = 0;
    int n_max_k_cuts = 0;
    int n_compat_cuts = 0;

    for (int snp = 0; snp < ref.L; ++snp) {
        PbwtStep step = state.advance(snp, block_start);
        bool target_ok = target_compat.constrain(ref, snp);
        bool over_max_k = step.n_patterns > max_k;
        if (!over_max_k && target_ok) {
            continue;
        }
        if (snp == block_start) {
            if (over_max_k) {
                throw std::invalid_argument("block_max_k is smaller than the number of alleles at SNP "
                                            + std::to_string(snp) + ".");
            }
            throw std::runtime_error("A target allele at reference SNP " + std::to_string(snp)
                                     + " is absent from the reference.");
        }

        snp_offsets.push_back(snp);
        n_max_k_cuts += over_max_k;
        n_compat_cuts += !target_ok;
        block_start = snp;
        target_compat.reset();
        if (!target_compat.constrain(ref, snp)) {
            throw std::runtime_error("A target allele at reference SNP " + std::to_string(snp)
                                     + " is absent from the reference.");
        }
        if (step.n_alleles > max_k) {
            throw std::invalid_argument("block_max_k is smaller than the number of alleles at SNP "
                                        + std::to_string(snp) + ".");
        }
    }
    snp_offsets.push_back(ref.L);
    return label_blocks(ref, std::move(snp_offsets), n_max_k_cuts, n_compat_cuts);
}
