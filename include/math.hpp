#pragma once
#include <optional>

double delta_Elogx(double mu, double sigma2, double a, double b, std::optional<double> Elogx = std::nullopt);

double delta_ElogGamma_x(double mu, double sigma2, double a, double b);

