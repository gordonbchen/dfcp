#include <cmath>
#include <stdexcept>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "math.hpp"
#include "elbo.hpp"
#include "math.hpp"
#include "util.hpp"


double betaln(double a, double b) {
    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

double calc_gammal_elbo_Ell(const HyperParams& HP, const Params& params, const Clusters& clusters, int l) {
    if (clusters.soft) {
        double Ell = clusters.rs[l].size() * (
            delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, 0.0)
            - HP.K * delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, 0.0)
        );
        for (Cluster* a : clusters.rs[l]) {
            Ell -= delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, a->n_obs);
            for (int k = 0; k < HP.K; ++k) {
                Ell += delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, a->nk[k]);
            }
        }
        return Ell;
    }

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

double normal_entropy(double sigma2) {
    if (sigma2 <= 0.0) { throw std::invalid_argument("sigma2 must be positive."); };
    constexpr double pi = 3.141592653589793238462643383279502884;
    return 0.5 * std::log(2.0 * pi * std::exp(1.0) * sigma2);
}

double calc_elbo(const HyperParams& HP, const Params& params, const Clusters& clusters) {
    // alpha.
    double elbo = delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0, 0.0);
    elbo -= delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0, HP.N);
    elbo += (clusters.nR + HP.tau_1) * params.mu_log_alpha;
    elbo -= HP.tau_2 * params.mu_alpha;
    elbo += HP.tau_1 * std::log(HP.tau_2) - std::lgamma(HP.tau_1);

    for (int l = 0; l < HP.L - 1; ++l) {
        // d.
        // -1s cancel out with d_l entropy.
        int nQl = clusters.qs[l].size();
        int nRl1 = clusters.rs[l].size() + clusters.rs[l+1].size();
        elbo += (nQl - nRl1 + HP.v_1) * params.mu_log_d[l];
        elbo += HP.v_2 * delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1.0, 1.0);
        elbo -= betaln(HP.v_1, HP.v_2);
        elbo -= nQl * delta_ElogGamma_x(params.mu_d[l], params.sigma2_d[l], -1.0, 1.0);

        for (auto* b : clusters.qs[l]) {
            elbo += delta_ElogGamma_x(params.mu_d[l], params.sigma2_d[l], -1.0, b->n);
        }

        // alpha and d term.
        double z = params.mu_alpha / params.mu_d[l];
        double dd2 = boost::math::digamma(z) * (2.0 * z / (params.mu_d[l] * params.mu_d[l]));
        dd2 += boost::math::trigamma(z) * (z * z / (params.mu_d[l] * params.mu_d[l]));
        elbo += delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0 / params.mu_d[l], 0.0);
        elbo += 0.5 * params.sigma2_d[l] * dd2;

        z += nQl;
        dd2 = boost::math::digamma(z) * (2.0 * params.mu_alpha / std::pow(params.mu_d[l], 3));
        dd2 += boost::math::trigamma(z) * (params.mu_alpha * params.mu_alpha / std::pow(params.mu_d[l], 4));
        elbo -= delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0 / params.mu_d[l], nQl);
        elbo -= 0.5 * params.sigma2_d[l] * dd2;
    }

    // gamma.
    elbo += HP.L * (HP.phi_1 * std::log(HP.phi_2) - std::lgamma(HP.phi_1));
    for (int l = 0; l < HP.L; ++l) {
        elbo += HP.phi_1*params.mu_log_gamma[l] - HP.phi_2*params.mu_gamma[l];
        elbo += calc_gammal_elbo_Ell(HP, params, clusters, l);
    }

    // Clusters.
    for (auto* a : clusters.rs[0]) {
        elbo += std::lgamma(static_cast<double>(a->n));
    }
    for (int l = 0; l < HP.L - 1; ++l) {
        for (auto* a : clusters.rs[l]) {
            elbo += std::lgamma(static_cast<double>(a->children.size()));
            elbo -= std::lgamma(static_cast<double>(a->n));
        }
        for (auto* a : clusters.rs[l + 1]) {
            elbo += std::lgamma(static_cast<double>(a->parents.size()));
        }
    }

    // Variational entropy.
    elbo += normal_entropy(params.sigma2_log_alpha);
    for (int l = 0; l < HP.L; ++l) {
        elbo += normal_entropy(params.sigma2_log_gamma[l]);
    }
    for (int l = 0; l < HP.L - 1; ++l) {
        elbo += normal_entropy(params.sigma2_logit_d[l]);
    }

    if (clusters.noisy) {
        // eps.
        double digamma_alpha = boost::math::digamma(params.alpha_eps);
        double digamma_beta = boost::math::digamma(params.beta_eps);
        double digamma_sum = boost::math::digamma(params.alpha_eps + params.beta_eps);
        double Elog_mismatch = digamma_alpha - digamma_sum;
        double Elog_match = digamma_beta - digamma_sum;
        elbo += (params.alpha_eps-1)*Elog_mismatch + (params.beta_eps-1)*Elog_match
            - (clusters.n_obs - clusters.n_matches) * std::log(HP.K-1.0)
            - betaln(HP.lambda_1, HP.lambda_2);

        // eps entropy.
        elbo += betaln(params.alpha_eps, params.beta_eps)
            - (params.alpha_eps-1)*digamma_alpha - (params.beta_eps-1)*digamma_beta
            + (params.alpha_eps+params.beta_eps-2) * digamma_sum;
    }
    return elbo;
}

