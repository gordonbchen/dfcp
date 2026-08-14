#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <cstdint>
#include <utility>


class EarlyStopping {
    private:
        double min_val = std::numeric_limits<double>::infinity();
        int steps_since_min = 0;

    public:
        int step = 0;
        const int patience;
        const bool minimize;
        const double tol;

        EarlyStopping(int patience_ = 3, bool minimize_ = true, double tol_ = 1e-5) :
            patience(patience_), minimize(minimize_), tol(tol_)
        {}

        void update(double x) {
            ++step;
            if (!minimize) {
                x = -x;
            }
            if (min_val-x > tol) {
                steps_since_min = 0;
                min_val = x;
                return;
            }
            ++steps_since_min;
        }

        bool converged() const {
            return steps_since_min >= patience;
        }
};

inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

std::vector<size_t> count_emissions(const std::vector<int8_t>& x, const int N, const int L, const int K);
std::vector<int8_t> get_emission_modes(const std::vector<size_t>& emission_counts, const int L, const int K);
std::pair<std::vector<int8_t>, std::vector<size_t>> count_minor_alleles(
    const std::vector<size_t>& emission_counts, const std::vector<int>& masked_ls, const int n_masked_ls, const int K
);

double parse_double(char *s);
int parse_int(char *s);

