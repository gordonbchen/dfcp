#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "clusters.hpp"
#include "hyperparams.hpp"
#include "model_array.hpp"
#include "params.hpp"
#include "seq_array.hpp"


double get_new_cluster_emission_ll(
    int emission, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
);


struct ScalarObs {
    const ModelArray* model;
    const SeqArray* target;
    const std::unordered_map<int, int>* obs_ls;
    int seq_idx;
    const Clusters& clusters;
    const Params& params;
    const HyperParams& HP;

    ScalarObs(
        const ModelArray& model_, int seq_idx_,
        const Clusters& clusters_, const Params& params_, const HyperParams& HP_
    );

    ScalarObs(
        const SeqArray& target_, int target_idx_,
        const std::unordered_map<int, int>& obs_ls_,
        const Clusters& clusters_, const Params& params_, const HyperParams& HP_
    );

    int emission(int l) const;

    double emission_ll(Cluster* a, int l) const {
        return a == nullptr ?
            get_new_cluster_emission_ll(emission(l), l, clusters, params, HP) : 0.0;
    }

    const std::vector<Cluster*>& matching_as(int l) const {
        return clusters.get_matching_as(l, emission(l));
    }

    bool matches(Cluster* a, int l) const {
        int observed_emission = emission(l);
        return observed_emission == -1 || a->emission == observed_emission;
    }
};


struct BlockObs {
    std::vector<std::vector<Cluster*>> compatible_clusters;
    std::vector<std::uint8_t> emission_is_compatible;
    std::vector<double> new_emission_ll;
    std::vector<double> new_cluster_obs_ll;
    std::vector<double> compatible_emission_ll;
    const HyperParams& HP;

    BlockObs(
        const ModelArray& ref, const SeqArray& target, int target_idx,
        const std::unordered_map<int, int>& obs_ls,
        const Clusters& clusters, const Params& params, const HyperParams& HP
    );

    double emission_ll(Cluster* a, int l) const {
        return a == nullptr ? new_cluster_obs_ll[l] : 0.0;
    }

    const std::vector<Cluster*>& matching_as(int l) const {
        return compatible_clusters[l];
    }

    bool matches(Cluster* a, int l) const {
        return is_emission_compatible(l, a->emission);
    }

    bool is_emission_compatible(int l, int emission) const {
        return emission_is_compatible[HP.emission_idx(l, emission)] != 0;
    }

    double new_emission_ll_at(int l, int emission) const {
        return new_emission_ll[HP.emission_idx(l, emission)];
    }
};
