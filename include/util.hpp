#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include "seq_array.hpp"


inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

inline int8_t get_xil(const SeqArray& x, int i, int l, const std::unordered_map<int, int> *obs_ls) {
    if (obs_ls == nullptr) {
        return x(i, l);
    }
    auto mapped_l = obs_ls->find(l);
    return mapped_l == obs_ls->end() ? -1 : x(i, mapped_l->second);
}

double parse_double(const char* value);
int parse_int(const char* value);

class EarlyStopping {
    double best;
    int steps_since_best = 0;

    public:
        int step = 0;
        const int patience;
        const double tol;
        const int max_steps;

        EarlyStopping(int patience, double tol, int max_steps) :
            best(-std::numeric_limits<double>::infinity()),
            patience(patience), tol(tol), max_steps(max_steps) {}

        void update(double value) {
            ++step;
            if (value - best > tol) {
                best = value;
                steps_since_best = 0;
                return;
            }
            ++steps_since_best;
        }

        bool converged() const {
            return step >= max_steps || steps_since_best >= patience;
        }
};
