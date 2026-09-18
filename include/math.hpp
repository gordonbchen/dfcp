#pragma once
#include <optional>
#include <vector>

double delta_Elogx(double mu, double sigma2, double a, double b, std::optional<double> Elogx = std::nullopt);

double delta_ElogGamma_x(double mu, double sigma2, double a, double b);

double log_sum_exp(double x1, double x2);

void normalize_ll(std::vector<double>& ll, int L, int K);
