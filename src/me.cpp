#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>
#include "hyperparams.hpp"
#include "impute_io.hpp"
#include "params.hpp"
#include "clusters.hpp"
#include "max.hpp"
#include "fwd_bkwd.hpp"
#include "expect.hpp"
#include "elbo.hpp"
#include "json.hpp"
#include "seq_array.hpp"
#include "util.hpp"


enum class InitMode { viterbi, block, pbwt };

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
    InitMode init_mode, int pbwt_match_len, bool pbwt_match_curr, bool init_only,
    const SeqArray& x_train,
    Json& json
) {
    // Init clusters.
    auto t0 = std::chrono::steady_clock::now();
    switch (init_mode) {
        case InitMode::block:
            clusters.block_init(x_train);
            break;
        case InitMode::pbwt:
            if (HP.K != 2) { throw std::invalid_argument("PBWT init only supported for K=2."); }
            clusters.pbwt_init(x_train, pbwt_match_len, pbwt_match_curr);
            break;
        case InitMode::viterbi:
            HP.N = 0;
            add_seqs(x_train, clusters, params, HP);
            break;
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
        expect_step(HP, params, clusters);
        auto t1 = std::chrono::steady_clock::now();
        max_step(x_train, clusters, params, HP);
        if (clusters.emit_mode == EmitMode::noisy) {
            max_cluster_emissions(clusters, params, HP);
        }
        auto t2 = std::chrono::steady_clock::now();
        elbo = calc_elbo(HP, params, clusters);
        auto t3 = std::chrono::steady_clock::now();
        early_stop.update(elbo);

        auto t_expect = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        auto t_max = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        auto t_elbo = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
        auto t_step = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t0).count();
        std::cerr << early_stop.step << ": elbo=" << elbo
            << " t_expect=" << t_expect << "ms t_max=" << t_max
            << "ms t_elbo=" << t_elbo << "ms t_step=" << t_step << "ms\n";
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
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    std::vector<Cluster*> viterbi_clusters{get_viterbi_clusters(x, i, &obs_ls, clusters, params, HP)};
    std::vector<double> viterbi_probs;
    int n_masked_ls = HP.L - obs_ls.size();
    viterbi_probs.reserve(n_masked_ls * HP.K);
    for (int l = 0; l < HP.L; ++l) {
        if (obs_ls.contains(l)) { continue; }
        for (int k = 0; k < HP.K; ++k) {
            double p = get_cluster_emission_ll(viterbi_clusters[2*l], k, l, clusters, params, HP);
            viterbi_probs.emplace_back(p);
        }
    }
    normalize_ll(viterbi_probs, n_masked_ls, HP.K);
    return viterbi_probs;
}

void impute(
    const SeqArray& x_val, const std::unordered_map<int, int>& obs_ls,
    bool viterbi,
    const Clusters& clusters, const Params& params, const HyperParams& HP,
    const char* prob_file, Json& json
) {
    int n_masked_ls = HP.L - static_cast<int>(obs_ls.size());
    ImputeProbWriter prob_writer(prob_file, x_val.N, n_masked_ls);

    std::chrono::steady_clock::duration impute_dur{};
    for (int i = 0; i < x_val.N; ++i) {
        auto t0 = std::chrono::steady_clock::now();
        std::vector<double> seq_probs = viterbi ?
            get_viterbi_impute_probs(x_val, i, obs_ls, clusters, params, HP)
            : fwd_bkwd(x_val, i, obs_ls, clusters, params, HP);
        prob_writer.write_row(seq_probs);
        impute_dur += std::chrono::steady_clock::now() - t0;
    }
    prob_writer.finish();

    auto t_impute = std::chrono::duration_cast<std::chrono::milliseconds>(impute_dur).count();
    std::cerr << "t_impute=" << t_impute << "ms\n";
    json.add("t_impute", t_impute);
}


std::unordered_map<int, int> read_obs_ls(const char *filename, int n_obs_ls, int n_total_loci) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open observed loci file.");
    }

    std::unordered_map<int, int> obs_ls;
    obs_ls.reserve(n_obs_ls);

    int i = 0;
    int previous_l = -1;
    int l;
    while (file >> l) {
        if (l < 0 || l >= n_total_loci) {
            throw std::runtime_error("Observed locus is outside the reference sequence.");
        }
        if (l <= previous_l) {
            throw std::runtime_error("Observed loci must be unique and strictly increasing.");
        }
        obs_ls.emplace(l, i);
        previous_l = l;
        ++i;
    }
    if (!file.eof()) {
        throw std::runtime_error("Observed loci file contains a non-integer value.");
    }
    if (i != n_obs_ls) {
        throw std::runtime_error("Observed loci file has the wrong number of rows.");
    }
    return obs_ls;
}

int main(int argc, char *argv[]) {
    Json json;

    // Read ref and target files.
    if (argc < 5) {
        throw std::invalid_argument("Requires ref bin, target obs bin, observed loci, and prob output file.");
    }
    std::cerr << "ref_file=" << argv[1] << " target_file=" << argv[2]
        << " observed_loci_file=" << argv[3] << " prob_file=" << argv[4] << '\n';
    json.add("ref_file", argv[1]).add("target_file", argv[2]).add("observed_loci_file", argv[3]);

    SeqArray x_train{read_seq_file(argv[1])};
    SeqArray x_val{read_seq_file(argv[2])};
    std::unordered_map<int, int> obs_ls{read_obs_ls(argv[3], x_val.L, x_train.L)};

    HyperParams HP{.N=x_train.N, .L=x_train.L, .K=2};

    // Parse optional args.
    EmitMode emit_mode = EmitMode::hard;

    InitMode init_mode = InitMode::pbwt;
    int pbwt_match_len = 5;
    bool pbwt_match_curr = true;
    bool init_only = false;

    bool viterbi_impute = false;

    int i = 5;
    while (i < argc) {
        if (i+1 >= argc) { throw std::invalid_argument("Arg has no value."); };

        std::string_view arg{argv[i]};
        if (arg == "--tau_1") { HP.tau_1 = parse_double(argv[i+1]); }
        else if (arg == "--tau_2") { HP.tau_2 = parse_double(argv[i+1]); }
        else if (arg == "--v_1") { HP.v_1 = parse_double(argv[i+1]); }
        else if (arg == "--v_2") { HP.v_2 = parse_double(argv[i+1]); }
        else if (arg == "--phi_1") { HP.phi_1 = parse_double(argv[i+1]); }
        else if (arg == "--phi_2") { HP.phi_2 = parse_double(argv[i+1]); }

        else if (arg == "--mode") {
            std::string_view value{argv[i+1]};
            if (value == "hard") { emit_mode = EmitMode::hard; }
            else if (value == "noisy") { emit_mode = EmitMode::noisy; }
            else if (value == "soft") { emit_mode = EmitMode::soft; }
            else { throw std::invalid_argument("mode must be hard, noisy, or soft."); }
        }
        else if (arg == "--lambda_1") { HP.lambda_1 = parse_double(argv[i+1]); }
        else if (arg == "--lambda_2") { HP.lambda_2 = parse_double(argv[i+1]); }

        else if (arg == "--init") {
            std::string_view value{argv[i+1]};
            if (value == "viterbi") { init_mode = InitMode::viterbi; }
            else if (value == "block") { init_mode = InitMode::block; }
            else if (value == "pbwt") { init_mode = InitMode::pbwt; }
            else { throw std::invalid_argument("init must be viterbi, block, or pbwt."); }
        }
        else if (arg == "--pbwt_match_len") { pbwt_match_len = parse_int(argv[i+1]); }
        else if (arg == "--pbwt_match_curr") { pbwt_match_curr = (parse_int(argv[i+1]) == 1); }
        else if (arg == "--init_only") { init_only = (parse_int(argv[i+1]) == 1); }

        else if (arg == "--viterbi_impute") { viterbi_impute = (parse_int(argv[i+1]) == 1); }

        else { throw std::invalid_argument("Arg not recognized."); }
        i += 2;
    }
    std::cerr << HP << '\n';

    // Init params and clusters.
    Params params{HP};
    Clusters clusters{HP, emit_mode};

    train_dfcp(
        clusters, params, HP,
        init_mode, pbwt_match_len, pbwt_match_curr, init_only,
        x_train,
        json
    );

    impute(
        x_val, obs_ls,
        viterbi_impute,
        clusters, params, HP,
        argv[4], json
    );
    std::cout << json.str() << '\n';
}
