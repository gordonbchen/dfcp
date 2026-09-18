#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "clusters.hpp"
#include "params.hpp"
#include "hyperparams.hpp"
#include "seq_array.hpp"
#include "math.hpp"
#include "obs.hpp"
#include "util.hpp"


double get_new_cluster_emission_ll(
    int emission, int l,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    if (emission == -1) {
        return 0.0;
    }

    int K = HP.n_emissions(l);
    int n_clusters_with_emission = clusters.rs_by_emit[HP.emission_idx(l, emission)].size();
    return delta_Elogx(
        params.mu_gamma[l], params.sigma2_gamma[l], 1.0,
        n_clusters_with_emission, params.mu_log_gamma[l]
    ) - delta_Elogx(
        params.mu_gamma[l], params.sigma2_gamma[l], K,
        clusters.rs[l].size(), params.mu_log_gamma[l]
    );
}



ScalarObs::ScalarObs(
    const ModelArray& model_, int seq_idx_,
    const Clusters& clusters_, const Params& params_, const HyperParams& HP_
) :
    model(&model_), target(nullptr), obs_ls(nullptr), seq_idx(seq_idx_),
    clusters(clusters_), params(params_), HP(HP_) {}

ScalarObs::ScalarObs(
    const SeqArray& target_, int target_idx_, const std::unordered_map<int, int>& obs_ls_,
    const Clusters& clusters_, const Params& params_, const HyperParams& HP_
) :
    model(nullptr), target(&target_), obs_ls(&obs_ls_), seq_idx(target_idx_),
    clusters(clusters_), params(params_), HP(HP_) {}

int ScalarObs::emission(int l) const {
    return model != nullptr ? static_cast<int>((*model)(seq_idx, l)) :
        get_xil(*target, seq_idx, l, obs_ls);
}


BlockObs::BlockObs(
    const ModelArray& ref, const SeqArray& target, int target_idx,
    const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) :
    compatible_clusters(HP.L),
    emission_is_compatible(HP.total_emissions(), 1),
    new_emission_ll(HP.total_emissions()),
    new_cluster_obs_ll(HP.L),
    compatible_emission_ll(HP.L),
    HP(HP)
{
    for (int l = 0; l < HP.L; ++l) {
        int K = HP.n_emissions(l);

        std::vector<std::pair<int, int>> observed_snps;
        for (int snp = ref.snp_start(l); snp < ref.snp_end(l); ++snp) {
            auto obs_it = obs_ls.find(snp);
            if (obs_it != obs_ls.end()) {
                observed_snps.emplace_back(snp, obs_it->second);
            }
        }

        for (int emission = 0; emission < K; ++emission) {
            std::size_t emission_idx = HP.emission_idx(l, emission);
            int representative = ref.rep(l, emission);
            for (auto [snp, target_snp] : observed_snps) {
                if (ref.snps(representative, snp) != target(target_idx, target_snp)) {
                    emission_is_compatible[emission_idx] = 0;
                    break;
                }
            }
            new_emission_ll[emission_idx] = get_new_cluster_emission_ll(
                emission, l, clusters, params, HP
            );
            if (emission_is_compatible[emission_idx] == 0) {
                continue;
            }
            const auto& emission_clusters = clusters.rs_by_emit[emission_idx];
            compatible_clusters[l].insert(
                compatible_clusters[l].end(), emission_clusters.begin(), emission_clusters.end()
            );
        }

        if (compatible_clusters[l].empty()) {
            throw std::runtime_error(
                "Target sequence " + std::to_string(target_idx) + ", block " + std::to_string(l)
                + " has no compatible reference pattern."
            );
        }

        double compatible_emissions_ll = -std::numeric_limits<double>::infinity();
        double all_emissions_ll = -std::numeric_limits<double>::infinity();
        for (int emission = 0; emission < K; ++emission) {
            std::size_t emission_idx = HP.emission_idx(l, emission);
            all_emissions_ll = log_sum_exp(all_emissions_ll, new_emission_ll[emission_idx]);
            if (emission_is_compatible[emission_idx] != 0) {
                compatible_emissions_ll = log_sum_exp(
                    compatible_emissions_ll, new_emission_ll[emission_idx]
                );
            }
        }
        compatible_emission_ll[l] = compatible_emissions_ll;
        new_cluster_obs_ll[l] = compatible_emissions_ll - all_emissions_ll;
    }
}
