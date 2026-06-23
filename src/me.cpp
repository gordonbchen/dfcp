#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <iomanip>
#include <string_view>
#include <random>
#include <cstdlib>
#include <stdexcept>
#include <chrono>
#include <tuple>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "max.hpp"
#include "expect.hpp"
#include "elbo.hpp"
#include "tree.hpp"
#include "util.hpp"


void parse_double(char *s, double& x) {
    char* end_ptr = nullptr;
    x = std::strtod(s, &end_ptr);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse arg value double."); };
}

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
            return steps_since_min > patience;
        }
};

struct SparseX {
    size_t i;
    size_t l;
    char x;
    SparseX(size_t i_, size_t l_, char x_) : i(i_), l(l_), x(x_) {}
};

int main(int argc, char *argv[]) {
    // Read seq file.
    if (argc < 2) { throw std::invalid_argument("Requires sequence file."); }
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "seq_file=" << argv[1] << '\n';
    std::ifstream seq_file(argv[1]);
    if (!seq_file.is_open()) { throw std::runtime_error("Failed to open sequence file."); };

    int N = 0;
    int L = 0;
    std::vector<char> x;
    std::string line;
    while (std::getline(seq_file, line)) {
        if (N == 0) {
            L = line.length();
        }
        for (char c : line) {
            if (c < '0' || c > '9') { throw std::runtime_error("Invalid allele char."); };
            x.push_back(c - '0');
        }
        ++N;
    }
    int K = *std::max_element(x.begin(), x.end()) + 1;
    HyperParams HP{.N=N, .L=L, .K=K};

    // Parse optional args.
    double val = -1.0;
    double mask = -1.0;
    char *tree_file_name = nullptr;

    int i = 2;
    while (i < argc) {
        if (i+1 >= argc) { throw std::invalid_argument("Arg has no value."); };

        std::string_view arg{argv[i]};
        if (arg == "--tau_1") { parse_double(argv[i+1], HP.tau_1); }
        else if (arg == "--tau_2") { parse_double(argv[i+1], HP.tau_2); }
        else if (arg == "--v_1") { parse_double(argv[i+1], HP.v_1); }
        else if (arg == "--v_2") { parse_double(argv[i+1], HP.v_2); }
        else if (arg == "--phi_1") { parse_double(argv[i+1], HP.phi_1); }
        else if (arg == "--phi_2") { parse_double(argv[i+1], HP.phi_2); }

        else if (arg == "--val") { parse_double(argv[i+1], val); }
        else if (arg == "--mask") { parse_double(argv[i+1], mask); }

        else if (arg == "--tree") { tree_file_name = argv[i+1]; }

        else { throw std::invalid_argument("Arg not recognized."); }
        i += 2;
    }

    std::cout << HP << '\n';
    std::cout << "val=" << val << ", mask=" << mask << '\n';

    // Split val for imputation.
    bool do_val = val > 0.0;
    if (do_val != (mask > 0.0)) { throw std::invalid_argument("If imputation val frac > 0, need mask frac > 0."); };
    int n_val_seqs = 0;
    std::vector<char> x_val_masked;
    x_val_masked.reserve(do_val ? static_cast<size_t>(val * x.size()) : 0);

    int n_masked_alleles = 0;
    std::vector<SparseX> x_val_true;
    x_val_true.reserve(do_val ? static_cast<size_t>(val * mask * x.size()) : 0);

    std::vector<char> x_train;
    x_train.reserve(static_cast<size_t>((1.0-val) * x.size()));

    if (val > 0.0) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::bernoulli_distribution val_dist(val);
        std::bernoulli_distribution mask_dist(mask);

        for (int i = 0; i < HP.N; ++i) {
            auto line = x.begin() + i*HP.L;
            if (!val_dist(gen)) {
                x_train.insert(x_train.end(), line, line + HP.L);
                continue;
            }
            x_val_masked.insert(x_val_masked.end(), line, line + HP.L);
            for (int l = 0; l < HP.L; ++l) {
                if (!mask_dist(gen)) { continue; }
                x_val_true.emplace_back(n_val_seqs, l, x_val_masked[idx2d(n_val_seqs, l, HP.L)]);
                x_val_masked[idx2d(n_val_seqs, l, HP.L)] = -1;
                ++n_masked_alleles;
            }
            ++n_val_seqs;
        }
        HP.N -= n_val_seqs;
        std::cout << HP << "\nn_val_seqs=" << n_val_seqs << " n_masked_alleles=" << n_masked_alleles << '\n';
    }
    else {
        x_train = std::move(x);
    }
    int n_train = HP.N;

    // Init params and clusters.
    Params params{HP};
    Clusters clusters{HP, x_train};

    EarlyStopping early_stop{2, false, 1e-3};
    double elbo = 0.0;
    while (!early_stop.converged()) {
        auto t0 = std::chrono::steady_clock::now();
        max_step(clusters, x_train, HP, params);
        auto t1 = std::chrono::steady_clock::now();
        expect_step(HP, params, clusters);
        auto t2 = std::chrono::steady_clock::now();
        elbo = calc_elbo(HP, params, clusters);
        auto t3 = std::chrono::steady_clock::now();
        early_stop.update(elbo);

        auto t_max = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        auto t_expect = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        auto t_elbo = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
        auto t_step = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t0).count();
        std::cout << early_stop.step << ": elbo=" << elbo
            << " t_max=" << t_max << "ms t_expect=" << t_expect << "ms t_elbo=" << t_elbo
            << "ms t_step=" << t_step << "ms\n";
    }

    // Impute.
    if (n_val_seqs > 0) {
        add_seqs(clusters, x_val_masked, params, HP);

        // TODO: modes should be of training data.
        std::vector<char> modes(HP.L);
        count_modes(modes, x_val_masked, n_val_seqs, HP.L, HP.K);

        int n_dfcp_correct = 0;
        int n_mode_correct = 0;
        for (SparseX& s : x_val_true) {
            if (s.x == clusters.r_assign[idx2d(n_train + s.i, s.l, HP.L)]->emission) {
                ++n_dfcp_correct;
            }
            if (s.x == modes[s.l]) {
                ++n_mode_correct;
            }
        }

        std::cout << "DFCP impute acc: " << n_dfcp_correct << '/' << n_masked_alleles << " = "
            << static_cast<double>(n_dfcp_correct) / static_cast<double>(n_masked_alleles) << '\n';
        std::cout << "Mode impute acc: " << n_mode_correct << '/' << n_masked_alleles << " = "
            << static_cast<double>(n_mode_correct) / static_cast<double>(n_masked_alleles) << '\n';
    }

    // Tree parsimony.
    if (tree_file_name != nullptr) {
        std::vector<std::unordered_map<int, std::tuple<int, int>>> coal_trees{parse_tree_file(tree_file_name, HP.L)};
        for (int l = 0; l < HP.L; ++l) {
            std::unordered_map<Cluster*, int> cluster_idxs;
            cluster_idxs.reserve(clusters.rs[l].size());
            int i = 0;
            for (Cluster* c : clusters.rs[l]) {
                cluster_idxs.emplace(c, i);
                ++i;
            }
            int excess_parsimony = calc_excess_parsimony(l, coal_trees[l], clusters, cluster_idxs);
            std::cout << "l=" << l << " excess_parsimony=" << excess_parsimony << '\n';
        }
    }
    return 0;
}

