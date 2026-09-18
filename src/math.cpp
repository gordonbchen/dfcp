#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>
#include <boost/math/special_functions/trigamma.hpp>
#include "math.hpp"
#include "util.hpp"


double delta_Elogx(double mu, double sigma2, double a, double b, std::optional<double> Elogx) {
    if ((Elogx.has_value()) && (b == 0.0)) {
        // TODO : negative a.
        return std::log(a) + Elogx.value();
    }

    double x = a*mu + b;
    return std::log(x) - 0.5*sigma2*a*a / (x*x);
}

double delta_ElogGamma_x(double mu, double sigma2, double a = 1.0, double b = 0.0) {
    double x = a*mu + b;
    return std::lgamma(x) + 0.5 * sigma2 * a * a * boost::math::trigamma(x);
}

double log_sum_exp(double x1, double x2) {
    if (x1 == -std::numeric_limits<double>::infinity()) {
        return x2;
    }
    if (x2 == -std::numeric_limits<double>::infinity()) {
        return x1;
    }
    double xmax = std::max(x1, x2);
    return xmax + std::log(std::exp(x1 - xmax) + std::exp(x2 - xmax));
}

void normalize_ll(std::vector<double>& ll, int L, int K) {
    for (int l = 0; l < L; ++l) {
        double sum_ll = ll[idx2d(l,0,K)];
        for (int k = 1; k < K; ++k) {
            sum_ll = log_sum_exp(sum_ll, ll[idx2d(l,k,K)]);
        }
        for (int k = 0; k < K; ++k) {
            ll[idx2d(l,k,K)] = std::exp(ll[idx2d(l,k,K)] - sum_ll);
        }
    }
}
