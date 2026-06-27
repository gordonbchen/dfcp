#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <limits>
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
#include "json.hpp"
#include "util.hpp"


void parse_double(char *s, double& x) {
    char* end_ptr = nullptr;
    x = std::strtod(s, &end_ptr);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse double arg value."); };
}

void parse_int(char *s, int& x) {
    char* end_ptr = nullptr;
    x = std::strtol(s, &end_ptr, 10);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse int arg value."); };
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
    Json json;

    // Read seq file.
    if (argc < 2) { throw std::invalid_argument("Requires sequence file."); }
    std::cerr << "seq_file=" << argv[1] << '\n';
    json.add("seq_file", argv[1]);
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

    char *tree_fname = nullptr;
    char *variant_pos_fname = nullptr;
    int variant_start_pos = -1;
    char *tree_vis_fname = nullptr;

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

        else if (arg == "--tree") { tree_fname = argv[i+1]; }
        else if (arg == "--variant_pos") { variant_pos_fname = argv[i+1]; }
        else if (arg == "--variant_start_pos") { parse_int(argv[i+1], variant_start_pos); }
        else if (arg == "--tree_vis") { tree_vis_fname = argv[i+1]; }

        else { throw std::invalid_argument("Arg not recognized."); }
        i += 2;
    }
    std::cerr << HP << '\n';

    // Split val for imputation.
    bool do_val = val > 0.0;
    if (do_val != (mask > 0.0)) { throw std::invalid_argument("If imputation val frac > 0, need mask frac > 0."); };
    if (do_val && (tree_fname != nullptr)) {throw std::invalid_argument("Does not support validation and tree eval.");}
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
        std::cerr << HP << "\nn_val_seqs=" << n_val_seqs << " n_masked_alleles=" << n_masked_alleles << '\n';
        json.add("n_val_seqs", n_val_seqs).add("n_masked_alleles", n_masked_alleles);
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

    std::vector<Json> train_log;
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
        std::cerr << early_stop.step << ": elbo=" << elbo
            << " t_max=" << t_max << "ms t_expect=" << t_expect << "ms t_elbo=" << t_elbo
            << "ms t_step=" << t_step << "ms\n";
        train_log.emplace_back();
        train_log[train_log.size()-1].add("elbo", elbo)
            .add("t_max", t_max).add("t_expect", t_expect).add("t_elbo", t_elbo).add("t_step", t_step);
    }
    json.add("train_log", train_log);

    // Impute.
    if (n_val_seqs > 0) {
        auto t0 = std::chrono::steady_clock::now();
        add_seqs(clusters, x_val_masked, params, HP);
        auto t1 = std::chrono::steady_clock::now();
        auto t_impute = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        std::vector<char> modes{count_modes(x_train, n_train, HP.L, HP.K)};

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
        double dfcp_impute_acc = static_cast<double>(n_dfcp_correct) / static_cast<double>(n_masked_alleles);
        double mode_impute_acc = static_cast<double>(n_mode_correct) / static_cast<double>(n_masked_alleles);

        std::cerr << "dfcp_impute_acc=" << dfcp_impute_acc << " t_impute=" << t_impute
            << "ms\nmode_impute_acc=" << mode_impute_acc << '\n';
        json.add("dfcp_impute_acc", dfcp_impute_acc).add("t_impute", t_impute).add("mode_impute_acc", mode_impute_acc);
    }

    // Tree parsimony.
    if (tree_fname != nullptr) {
        if ((variant_pos_fname == nullptr) || (variant_start_pos < 0)) {
            throw std::invalid_argument("Evaluating trees requires variant position file and variant start pos.");
        }
        std::vector<int> variant_pos{parse_pos_file_idx(variant_pos_fname, variant_start_pos, HP.L)};
        auto [coal_trees, recomb_pos] = parse_tree_file(tree_fname);
        std::vector<int> tree_idxs{get_tree_idxs(variant_pos, recomb_pos)};

        auto t0 = std::chrono::steady_clock::now();
        int excess_parsimony = 0;
        int emission_excess_parsimony = 0;
        for (int l = 0; l < HP.L; ++l) {
            std::unordered_map<Cluster*, int> cluster_idxs;
            cluster_idxs.reserve(clusters.rs[l].size());
            for (Cluster* c : clusters.rs[l]) {
                cluster_idxs.emplace(c, cluster_idxs.size());
            }
            std::vector<int> cluster_assignments(HP.N);
            for (int i = 0; i < HP.N; ++i) {
                cluster_assignments[i] = cluster_idxs.at(clusters.r_assign[idx2d(i, l, HP.L)]);
            }
            excess_parsimony += calc_excess_parsimony(
                coal_trees[tree_idxs[l]], cluster_assignments, clusters.rs[l].size()
            );

            std::vector<int> emission_clusters(HP.N);
            for (int i = 0; i < HP.N; ++i) {
                emission_clusters[i] = x_train[idx2d(i, l, HP.L)];
            }
            emission_excess_parsimony += calc_excess_parsimony(coal_trees[tree_idxs[l]], emission_clusters, -1);

            if ((tree_vis_fname != nullptr) && (l < 5)) {
                tree_to_dot(tree_vis_fname, coal_trees[tree_idxs[l]], cluster_assignments, emission_clusters, l, 5);
            }
        }
        auto t1 = std::chrono::steady_clock::now();
        auto t_parsimony = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

        double mean_excess_parsimony = static_cast<double>(excess_parsimony) / HP.L;
        double mean_emission_excess_parsimony = static_cast<double>(emission_excess_parsimony) / HP.L;
        std::cerr << "mean_excess_parsimony=" << mean_excess_parsimony
            << " t_parsimony=" << t_parsimony << "ms\n"
            << "mean_emission_excess_parsimony=" << mean_emission_excess_parsimony << '\n';
        json.add("mean_excess_parsimony", mean_excess_parsimony).add("t_parsimony", t_parsimony)
            .add("mean_emission_excess_parsimony", mean_emission_excess_parsimony);
    }
    std::cout << json.str() << '\n';
    return 0;
}

