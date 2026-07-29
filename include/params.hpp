#pragma once

#include <cmath>
#include <vector>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include "hyperparams.hpp"


struct Params {
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

    Params(const HyperParams& HP) :
        mu_alpha(HP.tau_1 / HP.tau_2),
        sigma2_alpha(HP.tau_1 / (HP.tau_2*HP.tau_2)),
        mu_log_alpha(boost::math::digamma(HP.tau_1) - std::log(HP.tau_2)),
        sigma2_log_alpha(boost::math::trigamma(HP.tau_1)),

        mu_d(HP.L-1, HP.v_1 / (HP.v_1 + HP.v_2)),
        sigma2_d(HP.L-1, (HP.v_1*HP.v_2) / (std::pow(HP.v_1+HP.v_2, 2) * (HP.v_1+HP.v_2+1))),
        mu_log_d(HP.L-1, boost::math::digamma(HP.v_1) - boost::math::digamma(HP.v_1 + HP.v_2)),
        sigma2_logit_d(HP.L-1, boost::math::trigamma(HP.v_1) + boost::math::trigamma(HP.v_2)),

        mu_gamma(HP.L, HP.phi_1 / HP.phi_2),
        sigma2_gamma(HP.L, HP.phi_1 / std::pow(HP.phi_2, 2)),
        mu_log_gamma(HP.L, boost::math::digamma(HP.phi_1) - std::log(HP.phi_2)),
        sigma2_log_gamma(HP.L, boost::math::trigamma(HP.phi_1)),

        alpha_eps(HP.lambda_1),
        beta_eps(HP.lambda_2)
    {}
};

