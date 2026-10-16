#pragma once

#include <cmath>
#include <vector>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include "hyperparams.hpp"


inline double calc_beta_mean(double alpha, double beta) {
    return alpha / (alpha + beta);
}

inline double calc_beta_var(double alpha, double beta) {
    return (alpha*beta) / (std::pow(alpha+beta, 2.0) * (alpha+beta+1.0));
}

inline double calc_gamma_mean(double shape, double rate) {
    return shape / rate;
}

inline double calc_gamma_var(double shape, double rate) {
    return shape / (rate * rate);
}


struct Params {
    int K;

    double mu_alpha;
    double sigma2_alpha;
    double mu_log_alpha;
    double sigma2_log_alpha;

    std::vector<double> mu_d;
    std::vector<double> sigma2_d;
    std::vector<double> mu_log_d;
    std::vector<double> sigma2_logit_d;

    std::vector<double> mu_gamma;
    std::vector<double> sigma2_gamma;
    std::vector<double> mu_log_gamma;
    std::vector<double> sigma2_log_gamma;

    double alpha_eps;
    double beta_eps;
    double mu_eps;
    double sigma2_eps;
    double Eeps_log_match;
    double Eeps_log_mismatch;

    Params(const HyperParams& HP) :
        K(HP.K),

        mu_alpha(calc_gamma_mean(HP.tau_1, HP.tau_2)),
        sigma2_alpha(calc_gamma_var(HP.tau_1, HP.tau_2)),
        mu_log_alpha(boost::math::digamma(HP.tau_1) - std::log(HP.tau_2)),
        sigma2_log_alpha(boost::math::trigamma(HP.tau_1)),

        mu_d(HP.L-1, calc_beta_mean(HP.v_1, HP.v_2)),
        sigma2_d(HP.L-1, calc_beta_var(HP.v_1, HP.v_2)),
        mu_log_d(HP.L-1, boost::math::digamma(HP.v_1) - boost::math::digamma(HP.v_1 + HP.v_2)),
        sigma2_logit_d(HP.L-1, boost::math::trigamma(HP.v_1) + boost::math::trigamma(HP.v_2)),

        mu_gamma(HP.L, calc_gamma_mean(HP.phi_1, HP.phi_2)),
        sigma2_gamma(HP.L, calc_gamma_var(HP.phi_1, HP.phi_2)),
        mu_log_gamma(HP.L, boost::math::digamma(HP.phi_1) - std::log(HP.phi_2)),
        sigma2_log_gamma(HP.L, boost::math::trigamma(HP.phi_1)),

        alpha_eps(HP.lambda_1),
        beta_eps(HP.lambda_2),
        mu_eps(calc_beta_mean(HP.lambda_1, HP.lambda_2)),
        sigma2_eps(calc_beta_var(HP.lambda_1, HP.lambda_2))
    {
        update_Eeps_log_match_mismatch();
    }

    void update_Eeps_log_match_mismatch() {
        double digamma_eps_sum = boost::math::digamma(alpha_eps + beta_eps);
        Eeps_log_match = boost::math::digamma(beta_eps) - digamma_eps_sum;
        Eeps_log_mismatch = boost::math::digamma(alpha_eps) - digamma_eps_sum - std::log(K-1.0);
    }
};

