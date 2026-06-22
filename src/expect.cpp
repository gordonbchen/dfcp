#include <cmath>
#include <stdexcept>
#include <omp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/logistic_sigmoid.hpp>
#include <boost/math/tools/minima.hpp>
#include "expect.hpp"
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "math.hpp"
#include "util.hpp"


template <typename F_nll, typename F_d2>
void laplace_log_approx(
    F_nll nll_log_func, F_d2 ll_log_d2_func,
    double& mu, double& sigma2,
    double& logx_mode, double& logx_var
) {
    // TODO: bits and max iter for find min.
    logx_mode = boost::math::tools::brent_find_minima(nll_log_func, -10.0, 10.0, 30).first;
    logx_var = -1.0 / ll_log_d2_func(std::exp(logx_mode));
    if (logx_var <= 0.0) { throw std::runtime_error("logx_var < 0."); };

    mu = std::exp(logx_mode + 0.5*logx_var);
    sigma2 = (std::exp(logx_var) - 1.0) * std::exp(2.0*logx_mode + logx_var);
    if (sigma2 <= 0.0) { throw std::runtime_error("sigma2 < 0."); };
}

double delta_ElogGamma_invx(double mu, double sigma2, double a, double b) {
    double x = a/mu + b;
    double d2 = boost::math::trigamma(x) * (a*a) / std::pow(mu, 4);
    d2 += boost::math::digamma(x) * 2*a / std::pow(mu, 3);
    return std::lgamma(x) + 0.5*sigma2 * d2;
}

double delta_ElogGamma_invx_d2_x(double mu, double sigma2, double a, double b) {
    double x = a/mu + b;
    double d2 = boost::math::trigamma(x) / (mu*mu);
    d2 += 0.5*sigma2 * (
        boost::math::polygamma(3, x) * (a*a) / std::pow(mu, 6)
        + boost::math::polygamma(2, x) * (6*a) / std::pow(mu, 5)
        + boost::math::polygamma(1, x) * 6 / std::pow(mu, 4)
    );
    return d2;
}

double ll_log_alpha(double log_alpha, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    double alpha = std::exp(log_alpha);

    double ll = std::lgamma(alpha) - std::lgamma(alpha + HP.N);
    ll += (clusters.nR + HP.tau_1 - 1)*std::log(alpha) - HP.tau_2*alpha;
    for (int l = 0; l < HP.L-1; ++l) {
        ll += delta_ElogGamma_invx(params.mu_d[l], params.sigma2_d[l], alpha, 0.0);
        ll -= delta_ElogGamma_invx(params.mu_d[l], params.sigma2_d[l], alpha, clusters.qs[l].size());
    }

    return ll + log_alpha;
}

double ll_log_alpha_d2(double alpha, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    double d2 = boost::math::trigamma(alpha) - boost::math::trigamma(alpha + HP.N);
    d2 += (1 - HP.tau_1 - clusters.nR) / (alpha*alpha);
    for (int l = 0; l < HP.L-1; ++l) {
        d2 += delta_ElogGamma_invx_d2_x(params.mu_d[l], params.sigma2_d[l], alpha, 0.0);
        d2 -= delta_ElogGamma_invx_d2_x(params.mu_d[l], params.sigma2_d[l], alpha, clusters.qs[l].size());
    }
    return alpha*alpha * d2 - 1;
}

double ll_log_gammal(double log_gamma, int l, const HyperParams& HP, const Clusters& clusters) {
    double gamma = std::exp(log_gamma);

    double ll = (HP.phi_1-1)*std::log(gamma) - HP.phi_2*gamma;
    ll += std::lgamma(HP.K*gamma) - std::lgamma(HP.K*gamma + clusters.rs[l].size());
    for (int k = 0; k < HP.K; ++k) {
        ll += std::lgamma(gamma + clusters.rs_by_emit[idx2d(l, k, HP.K)].size());
    }
    ll -= HP.K * std::lgamma(gamma);

    return ll + log_gamma;
}

double ll_log_gammal_d2(double gamma, int l, const HyperParams& HP, const Clusters& clusters) {
    double d2 = (1.0 - HP.phi_1) / (gamma*gamma);
    d2 += HP.K*HP.K * (boost::math::trigamma(HP.K*gamma) - boost::math::trigamma(HP.K*gamma + clusters.rs[l].size()));
    for (int k = 0; k < HP.K; ++k) {
        d2 += boost::math::trigamma(gamma + clusters.rs_by_emit[idx2d(l, k, HP.K)].size());
        d2 -= boost::math::trigamma(gamma);
    }
    return gamma*gamma * d2 - 1;
}

double delta_ElogGamma_x_d2_invx(double mu, double sigma2, double a, double b) {
    double x = a*mu + b;
    double d2 = boost::math::trigamma(x) * mu*mu * std::pow(a, 4);
    d2 += boost::math::digamma(x) * 2*mu * std::pow(a, 3);
    d2 += 0.5 * sigma2 * (
        boost::math::polygamma(3, x) * mu*mu * std::pow(a, 6)
        + boost::math::polygamma(2, x) * 6*mu * std::pow(a, 5)
        + boost::math::trigamma(x) * 6 * std::pow(a, 4)
    );
    return d2;
}

double ll_logit_dl(double logit_dl, int l, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    double d = boost::math::logistic_sigmoid(logit_dl);

    int nQl = clusters.qs[l].size();
    int nRl1 = clusters.rs[l].size() + clusters.rs[l+1].size();
    double ll = (nQl - nRl1 + HP.v_1 - 1.0) * std::log(d);
    ll += (HP.v_2 - 1.0) * std::log(1.0 - d);
    ll -= nQl * std::lgamma(1.0 - d);
    for (Cluster* b : clusters.qs[l]) {
        ll += std::lgamma(b->n - d);
    }
    ll += delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0/d, 0.0);
    ll -= delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0/d, nQl);

    return ll + std::log(d) + std::log(1.0-d);
}

double ll_logit_dl_d2(double d, int l, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    int nQl = clusters.qs[l].size();
    int nRl1 = clusters.rs[l].size() + clusters.rs[l+1].size();
    double d2 = (nRl1 - nQl + 1.0 - HP.v_1) / (d*d);
    d2 -= (HP.v_2 - 1.0) / std::pow(1.0 - d, 2);
    d2 -= nQl * boost::math::trigamma(1.0 - d);
    for (Cluster* b : clusters.qs[l]) {
        d2 += boost::math::trigamma(b->n - d);
    }

    d2 += delta_ElogGamma_x_d2_invx(params.mu_alpha, params.sigma2_alpha, 1.0/d, 0.0);
    d2 -= delta_ElogGamma_x_d2_invx(params.mu_alpha, params.sigma2_alpha, 1.0/d, nQl);

    d2 *= std::pow(d * (1.0-d), 2);
    d2 -= std::pow(1.0 - 2.0*d, 2) + 2.0*d*(1.0-d);
    return d2;
}

void expect_step(const HyperParams& HP, Params& params, const Clusters& clusters) {
    // alpha update.
    laplace_log_approx(
        [&](double log_alpha) -> double { return -ll_log_alpha(log_alpha, HP, params, clusters); },
        [&](double alpha) -> double { return ll_log_alpha_d2(alpha, HP, params, clusters); },
        params.mu_alpha, params.sigma2_alpha,
        params.mu_log_alpha, params.sigma2_log_alpha
    );

    #pragma omp parallel for
    for (int l = 0; l < HP.L; ++l) {
        // gamma_l update.
        laplace_log_approx(
            [&](double log_gammal) { return -ll_log_gammal(log_gammal, l, HP, clusters); },
            [&](double gammal) { return ll_log_gammal_d2(gammal, l, HP, clusters); },
            params.mu_gamma[l], params.sigma2_gamma[l],
            params.mu_log_gamma[l], params.sigma2_log_gamma[l]
        );
    }

    #pragma omp parallel for
    for (int l = 0; l < HP.L-1; ++l) {
        // d_l update in logit space.
        double logit_dl_mode = boost::math::tools::brent_find_minima(
            [&](double logit_dl) {return -ll_logit_dl(logit_dl, l, HP, params, clusters); },
            -10.0, 10.0, 30
        ).first;
        double dl_mode = boost::math::logistic_sigmoid(logit_dl_mode);
        params.sigma2_logit_d[l] = -1.0 / ll_logit_dl_d2(dl_mode, l, HP, params, clusters);
        if (params.sigma2_logit_d[l] <= 0.0) { throw std::runtime_error("sigma2_logit_d < 0."); };

        params.mu_d[l] = dl_mode + 0.5*params.sigma2_logit_d[l] * (1.0 - 2.0*dl_mode)*dl_mode*(1.0 - dl_mode);
        params.sigma2_d[l] = (
            dl_mode*dl_mode
            + 0.5*params.sigma2_logit_d[l] * (4.0*dl_mode - 6.0*dl_mode*dl_mode) * dl_mode * (1.0 - dl_mode)
            - params.mu_d[l]*params.mu_d[l]
        );
        if (params.sigma2_d[l] <= 0.0) { throw std::runtime_error("sigma2_d < 0."); };
        params.mu_log_d[l] = delta_Elogx(params.mu_d[l], params.sigma2_d[l], 1.0, 0.0);
    }
}

