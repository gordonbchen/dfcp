#include <vector>
#include <cstdint>
#include <algorithm>
#include <utility>
#include "pbwt.hpp"
#include "util.hpp"


std::pair<std::vector<int>, std::vector<int>> pbwt(const std::vector<int8_t>& x, int N, int L) {
    std::vector<int> a(N*(L+1), 0);
    std::vector<int> d(N*(L+1), 0);
    for (int i = 0; i < N; ++i) {
        a[idx2d(i,0,L+1)] = i;
    }

    int zeros_idx = 0;
    int ones_idx = 0;
    std::vector<int> zeros(N, 0);
    std::vector<int> ones(N, 0);
    std::vector<int> zeros_d(N, 0);
    std::vector<int> ones_d(N, 0);

    for (int l = 0; l < L; ++l) {
        zeros_idx = 0;
        ones_idx = 0;

        int p = l + 1;
        int q = l + 1;

        for (int i = 0; i < N; ++i) {
            p = std::max(p, d[idx2d(i,l,L+1)]);
            q = std::max(q, d[idx2d(i,l,L+1)]);

            int ai = a[idx2d(i,l,L+1)];
            if (x[idx2d(ai,l,L)] == 0) {
                zeros[zeros_idx] = ai;
                zeros_d[zeros_idx] = p;
                p = 0;
                ++zeros_idx;
            }
            else {
                ones[ones_idx] = ai;
                ones_d[ones_idx] = q;
                q = 0;
                ++ones_idx;
            }
        }

        for (int i = 0; i < zeros_idx; ++i) {
            a[idx2d(i,l+1,L+1)] = zeros[i];
            d[idx2d(i,l+1,L+1)] = zeros_d[i];
        }
        for (int i = 0; i < ones_idx; ++i) {
            a[idx2d(zeros_idx+i,l+1,L+1)] = ones[i];
            d[idx2d(zeros_idx+i,l+1,L+1)] = ones_d[i];
        }
    }

    return {std::move(a), std::move(d)};
}

