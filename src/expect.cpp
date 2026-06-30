#include <cmath>
#include <omp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include "expect.hpp"
#include "hyperparams.hpp"
#include "clusters.hpp"
#include "soft_clusters.hpp"
#include "util.hpp"


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


double ll_log_gammal(double log_gamma, int l, const HyperParams& HP, const SoftClusters& clusters) {
    double gamma = std::exp(log_gamma);

    double ll = (HP.phi_1-1)*std::log(gamma) - HP.phi_2*gamma;
    ll += clusters.rs[l].size() * (std::lgamma(HP.K*gamma) + HP.K*std::lgamma(gamma));

    for (SoftCluster* a : clusters.rs[l]) {
        ll -= std::lgamma(HP.K*gamma + a->n);
        for (int k = 0; k < HP.K; ++k) {
            ll += std::lgamma(gamma + a->nk[k]);
        }
    }
    return ll + log_gamma;
}

double ll_log_gammal_d2(double gamma, int l, const HyperParams& HP, const SoftClusters& clusters) {
    double d2 = (1.0 - HP.phi_1) / (gamma*gamma);
    d2 += clusters.rs[l].size() * (
        HP.K*HP.K * boost::math::trigamma(HP.K*gamma) - HP.K*boost::math::trigamma(gamma)
    );

    for (SoftCluster* a : clusters.rs[l]) {
        d2 -= HP.K*HP.K * boost::math::trigamma(HP.K*gamma + a->n);
        for (int k = 0; k < HP.K; ++k) {
            d2 += boost::math::trigamma(gamma + a->nk[k]);
        }
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

