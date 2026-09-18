#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>


struct HyperParams {
    int N;
    int L;
    std::vector<int> n_emissions_by_locus;
    std::vector<std::size_t> emission_offsets;
    double tau_1 = 1.0;
    double tau_2 = 1.0;
    double v_1 = 1.0;
    double v_2 = 1.0;
    double phi_1 = 2.0;
    double phi_2 = 2.0;

    HyperParams(int N_, std::vector<int> n_emissions) :
        N(N_),
        L(static_cast<int>(n_emissions.size())),
        n_emissions_by_locus(std::move(n_emissions)),
        emission_offsets(static_cast<std::size_t>(L) + 1, 0)
    {
        for (int l = 0; l < L; ++l) {
            emission_offsets[l + 1] = emission_offsets[l] + n_emissions_by_locus[l];
        }
    }

    int n_emissions(int l) const {
        return n_emissions_by_locus[l];
    }

    std::size_t emission_idx(int l, int emission) const {
        return emission_offsets[l] + emission;
    }

    std::size_t total_emissions() const {
        return emission_offsets.back();
    }
};

inline std::ostream& operator<<(std::ostream& os, const HyperParams& HP) {
    auto [min_k, max_k] = std::minmax_element(
        HP.n_emissions_by_locus.begin(), HP.n_emissions_by_locus.end()
    );
    return (
        os << "HP: N=" << HP.N << ", L=" << HP.L
        << ", K_min=" << *min_k << ", K_max=" << *max_k
        << ", tau_1=" << HP.tau_1 << ", tau_2=" << HP.tau_2
        << ", v_1=" << HP.v_1 << ", v_2=" << HP.v_2
        << ", phi_1=" << HP.phi_1 << ", phi_2=" << HP.phi_2
    );
}
