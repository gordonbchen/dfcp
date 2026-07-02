#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <cstdint>


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

struct SparseX {
    size_t i;
    size_t l;
    int8_t x;
    SparseX(size_t i_, size_t l_, int8_t x_) : i(i_), l(l_), x(x_) {}
};

inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

std::vector<int8_t> count_modes(const std::vector<int8_t>& x, const int N, const int L, const int K);

double parse_double(char *s);
int parse_int(char *s);

