#include <cmath>
#include <boost/math/special_functions/trigamma.hpp>
#include <optional>
#include "math.hpp"


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

