#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "soft_clusters.hpp"
#include "elbo.hpp"
#include "math.hpp"
#include "util.hpp"


double calc_gammal_elbo_Ell(const HyperParams& HP, const Params& params, const Clusters& clusters, int l) {
    double Ell = delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, 0.0);
    Ell -= delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, clusters.rs[l].size());
    for (int k = 0; k < HP.K; ++k) {
        Ell += delta_ElogGamma_x(
            params.mu_gamma[l], params.sigma2_gamma[l], 1.0, clusters.rs_by_emit[idx2d(l, k, HP.K)].size()
        );
    }
    Ell -= HP.K * delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, 0.0);
    return Ell;
}

double calc_gammal_elbo_Ell(const HyperParams& HP, const Params& params, const SoftClusters& clusters, int l) {
    double Ell = clusters.rs[l].size() * (
        delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, 0.0)
        - HP.K * delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, 0.0)
    );
    for (SoftCluster* a : clusters.rs[l]) {
        Ell -= delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n);
        for (int k = 0; k < HP.K; ++k) {
            Ell += delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, a->nk[k]);
        }
    }
    return Ell;
}

