#include <algorithm>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <string_view>
#include <cstdlib>
#include <stdexcept>
#include <chrono>
#include <limits>
#include <unordered_map>
#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "max.hpp"
#include "fwd_bkwd.hpp"
#include "expect.hpp"
#include "elbo.hpp"
#include "json.hpp"


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

void train_dfcp(
    Clusters& clusters, Params& params, HyperParams& HP,

    bool block_init,
    bool pbwt_init, int pbwt_match_len, bool pbwt_match_curr,
    bool init_only,

    const std::vector<int8_t>& x_train,

    Json& json
) {
    // Init clusters.
    auto t0 = std::chrono::steady_clock::now();
    if (block_init) {
        clusters.block_init(x_train);
    }
    else if (pbwt_init) {
        if (HP.K != 2) { throw std::invalid_argument("pbwt_init only supported for K=2."); }
        clusters.pbwt_init(x_train, pbwt_match_len, pbwt_match_curr);
    }
    else {
        int N = HP.N;
        HP.N = 0;
        add_seqs(x_train.begin(), N, clusters, params, HP);
    }
    auto t1 = std::chrono::steady_clock::now();
    auto t_init = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cerr << "t_init=" << t_init << '\n';
    json.add("t_init", t_init);

    if (init_only) { return; }

    // Train.
    EarlyStopping early_stop{2, false, 1e-3};
    double elbo = 0.0;

    std::vector<Json> train_log;
    while (!early_stop.converged()) {
        auto t0 = std::chrono::steady_clock::now();
        max_step(x_train, clusters, params, HP);
        if (clusters.noisy) {
            max_cluster_emissions(clusters, params, HP);
        }
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

    Json param_log;
    param_log.add("mu_alpha", params.mu_alpha).add("mu_gamma", params.mu_gamma).add("mu_d", params.mu_d)
        .add("alpha_eps", params.alpha_eps).add("beta_eps", params.beta_eps);
    json.add("params", param_log);
}


std::vector<double> get_viterbi_impute_probs(
    std::vector<int8_t>::const_iterator xi, const std::unordered_map<int, int>& unmasked_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<Cluster*> viterbi_clusters{get_viterbi_clusters(xi, &unmasked_ls, clusters, params, HP)};
    std::vector<double> viterbi_probs;
    int n_masked_ls = HP.L - unmasked_ls.size();
    viterbi_probs.reserve(n_masked_ls * HP.K);
    for (int l = 0; l < HP.L; ++l) {
        if (unmasked_ls.contains(l)) { continue; }
        for (int k = 0; k < HP.K; ++k) {
            double p = get_cluster_emission_ll(viterbi_clusters[2*l], k, l, clusters, params, HP);
            viterbi_probs.emplace_back(p);
        }
    }
    normalize_ll(viterbi_probs, n_masked_ls, HP.K);
    return viterbi_probs;
}

void impute(
    const std::vector<int8_t>& x_val, int n_val_seqs, const std::unordered_map<int, int>& unmasked_ls,
    bool viterbi,
    const Clusters& clusters, const Params& params, const HyperParams& HP,
    Json& json
) {
    Json prob_log;
    std::vector<std::vector<double>> probs;
    probs.reserve(n_val_seqs);

    std::chrono::steady_clock::duration impute_dur{};
    for (int i = 0; i < n_val_seqs; ++i) {
        auto xi = x_val.begin() + i*HP.L;
        auto t0 = std::chrono::steady_clock::now();
        std::vector<double> seq_probs = viterbi ?
            get_viterbi_impute_probs(xi, unmasked_ls, clusters, params, HP)
            : fwd_bkwd(xi, unmasked_ls, clusters, params, HP);
        impute_dur += std::chrono::steady_clock::now() - t0;

        probs.emplace_back(std::move(seq_probs));
    }
    json.add("probs", probs);

    auto t_impute = std::chrono::duration_cast<std::chrono::milliseconds>(impute_dur).count();
    std::cerr << "t_impute=" << t_impute << "ms\n";
    json.add("t_impute", t_impute);
}


std::vector<int8_t> read_seq_file(char *filename, int& N, int& L, int& K) {
    std::ifstream file(filename);
    if (!file.is_open()) { throw std::runtime_error("Failed to open seq file."); };

    std::vector<int8_t> x;

    std::string line;
    N = 0;
    while (std::getline(file, line)) {
        if (N == 0) {
            L = line.size();
        }

        for (int l = 0; l < L; ++l) {
            if ((line[l] >= '0') && (line[l] <= '9')) {
                x.push_back(line[l] - '0');
                K = std::max(K, (line[l] - '0') + 1);
            }
            else {
                throw std::runtime_error("Invalid allele char.");
            };
        }
        ++N;
    }
    return x;
}

std::unordered_map<int, int> read_unmasked_ls(char *filename, int n_unmasked_ls) {
    std::ifstream file(filename);
    if (!file.is_open()) { throw std::runtime_error("Failed to open unmasked ls file."); };

    std::unordered_map<int, int> unmasked_ls;
    unmasked_ls.reserve(n_unmasked_ls);

    int i = 0;
    int l;
    while (file >> l) {
        unmasked_ls[l] = i;
        ++i;
        if (i > n_unmasked_ls) { break; }
    }
    if (unmasked_ls.size() != static_cast<size_t>(n_unmasked_ls)) {
        throw std::runtime_error("did not read n_unmasked_ls from file.");
    }
    return unmasked_ls;
}

double parse_double(char *s) {
    char* end_ptr = nullptr;
    double x = std::strtod(s, &end_ptr);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse double arg value."); };
    return x;
}

int parse_int(char *s) {
    char* end_ptr = nullptr;
    int x = std::strtol(s, &end_ptr, 10);
    if (end_ptr == s) { throw std::invalid_argument("Failed to parse int arg value."); };
    return x;
}

int main(int argc, char *argv[]) {
    Json json;

    // Read ref and target files.
    if (argc < 4) { throw std::invalid_argument("Requires ref and target seq files, and masked ls file."); }
    std::cerr << "ref_file=" << argv[1] << " target_file=" << argv[2]
        << " unmasked_ls_file=" << argv[3] << '\n';
    json.add("ref_file", argv[1]).add("target_file", argv[2]).add("unmasked_ls_file", argv[3]);

    int n_train_seqs;
    int L;
    int K = 2;
    std::vector<int8_t> x_train{read_seq_file(argv[1], n_train_seqs, L, K)};

    int n_val_seqs;
    int n_unmasked_ls;
    std::vector<int8_t> x_val{read_seq_file(argv[2], n_val_seqs, n_unmasked_ls, K)};

    std::unordered_map<int, int> unmasked_ls{read_unmasked_ls(argv[3], n_unmasked_ls)};

    HyperParams HP{.N=n_train_seqs, .L=L, .K=K};

    // Parse optional args.
    bool noisy = false;
    bool soft = false;

    bool block_init = false;
    bool pbwt_init = false;
    int pbwt_match_len = 5;
    bool pbwt_match_curr = true;
    bool init_only = false;

    bool viterbi_impute = false;

    int i = 2;
    while (i < argc) {
        if (i+1 >= argc) { throw std::invalid_argument("Arg has no value."); };

        std::string_view arg{argv[i]};
        if (arg == "--tau_1") { HP.tau_1 = parse_double(argv[i+1]); }
        else if (arg == "--tau_2") { HP.tau_2 = parse_double(argv[i+1]); }
        else if (arg == "--v_1") { HP.v_1 = parse_double(argv[i+1]); }
        else if (arg == "--v_2") { HP.v_2 = parse_double(argv[i+1]); }
        else if (arg == "--phi_1") { HP.phi_1 = parse_double(argv[i+1]); }
        else if (arg == "--phi_2") { HP.phi_2 = parse_double(argv[i+1]); }

        else if (arg == "--noisy") { noisy = parse_int(argv[i+1]) == 1; }
        else if (arg == "--lambda_1") { HP.lambda_1 = parse_double(argv[i+1]); }
        else if (arg == "--lambda_2") { HP.lambda_2 = parse_double(argv[i+1]); }

        else if (arg == "--soft") { soft = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--block_init") { block_init = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--pbwt_init") { pbwt_init = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--pbwt_match_len") { pbwt_match_len = parse_int(argv[i+1]); }
        else if (arg == "--pbwt_match_curr") { pbwt_match_curr = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--init_only") { init_only = (parse_int(argv[i+1]) == 1); }

        else if (arg == "--viterbi_impute") { viterbi_impute = (parse_int(argv[i+1]) == 1); }

        else { throw std::invalid_argument("Arg not recognized."); }
        i += 2;
    }
    if (noisy && soft) { throw std::invalid_argument("noisy cannot be combined with soft."); }
    if (block_init && pbwt_init) { throw std::invalid_argument("cannot do block and pbwt init."); }
    std::cerr << HP << '\n';

    // Init params and clusters.
    Params params{HP};
    Clusters clusters{HP, soft, noisy};

    train_dfcp(
        clusters, params, HP,

        block_init,
        pbwt_init, pbwt_match_len, pbwt_match_curr,
        init_only,

        x_train,

        json
    );

    impute(
        x_val, n_val_seqs, unmasked_ls,
        viterbi_impute,
        clusters, params, HP,
        json
    );
    std::cout << json.str() << '\n';
}

