#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "impute_io.hpp"
#include "seq_array.hpp"


int main(int argc, char* argv[]) {
    if (argc != 4) {
        throw std::invalid_argument("Usage: eval_impute PROBS.bin MASKED_TRUTH.bin OUTPUT.bin");
    }

    ImputeProbReader prob_reader(argv[1]);
    const ImputeProbHeader& header = prob_reader.header();

    SeqArray x_val_true{read_seq_file(argv[2])};
    if (header.n_sequences != static_cast<std::uint32_t>(x_val_true.N)
        || header.n_loci != static_cast<std::uint32_t>(x_val_true.L)) {
        throw std::runtime_error("Probability and truth dimensions do not match.");
    }

    std::size_t n_loci = header.n_loci;
    std::vector<double> sum_q(n_loci, 0.0);
    std::vector<double> sum_q2(n_loci, 0.0);
    std::vector<double> sum_qy(n_loci, 0.0);
    std::vector<std::uint32_t> sum_y(n_loci, 0);
    std::vector<std::uint32_t> n_correct(n_loci, 0);
    std::vector<std::uint16_t> prob_row(n_loci);

    for (int i = 0; i < x_val_true.N; ++i) {
        prob_reader.read_row(prob_row);
        for (int l = 0; l < x_val_true.L; ++l) {
            std::uint16_t q = prob_row[l];
            int y = x_val_true(i, l);
            sum_q[l] += q;
            sum_q2[l] += static_cast<double>(q) * q;
            sum_qy[l] += static_cast<double>(q) * y;
            sum_y[l] += y;
            n_correct[l] += (static_cast<int>(q >= std::uint16_t{1} << 15) == y);
        }
    }
    std::vector<float> r2(n_loci, -1.0F);
    std::vector<float> accuracy(n_loci);
    double r2_sum = 0.0;
    double accuracy_sum = 0.0;
    std::size_t n_defined_r2 = 0;
    double n = x_val_true.N;
    for (std::size_t l = 0; l < n_loci; ++l) {
        double centered_q = n * sum_q2[l] - sum_q[l] * sum_q[l];
        double centered_y = n * sum_y[l] - static_cast<double>(sum_y[l]) * sum_y[l];
        double covariance = n * sum_qy[l] - sum_q[l] * sum_y[l];
        if (centered_q > 0.0 && centered_y > 0.0) {
            r2[l] = static_cast<float>(std::clamp(
                covariance * covariance / (centered_q * centered_y), 0.0, 1.0
            ));
            r2_sum += r2[l];
            ++n_defined_r2;
        }
        accuracy[l] = static_cast<float>(n_correct[l] / n);
        accuracy_sum += accuracy[l];
    }

    double mean_r2 = n_defined_r2 ? r2_sum / n_defined_r2 : -1.0;
    double mean_accuracy = accuracy_sum / n_loci;
    std::cerr << "n_sequences=" << x_val_true.N << " n_loci=" << x_val_true.L
        << " mean_r2=" << mean_r2 << " mean_accuracy=" << mean_accuracy
        << " undefined_r2=" << n_loci - n_defined_r2 << '\n';
    write_impute_eval_file(argv[3], r2, accuracy);
    return 0;
}
