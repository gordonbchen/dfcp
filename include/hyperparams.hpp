#pragma once

#include <iostream>


struct HyperParams {
    int N;
    int L;
    int K;
    double tau_1 = 1.0;
    double tau_2 = 1.0;
    double v_1 = 1.0;
    double v_2 = 1.0;
    double phi_1 = 2.0;
    double phi_2 = 2.0;
    double lambda_1 = 0.1;
    double lambda_2 = 1.0;
};

inline std::ostream& operator<<(std::ostream& os, const HyperParams& HP) {
    return (
        os << "HP: N=" << HP.N << ", L=" << HP.L << ", K=" << HP.K
        << ", tau_1=" << HP.tau_1 << ", tau_2=" << HP.tau_2
        << ", v_1=" << HP.v_1 << ", v_2=" << HP.v_2
        << ", phi_1=" << HP.phi_1 << ", phi_2=" << HP.phi_2
        << ", lambda_1=" << HP.lambda_1 << ", lambda_2=" << HP.lambda_2
    );
}

