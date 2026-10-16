#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <CLI/CLI.hpp>
#include "hyperparams.hpp"
#include "r_assign_io.hpp"
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
        const int max_steps;

        EarlyStopping(int patience_, bool minimize_, double tol_, int max_steps_) :
            patience(patience_), minimize(minimize_), tol(tol_), max_steps(max_steps_)
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
            return (step >= max_steps) || (steps_since_min >= patience);
        }
};

void train_dfcp(
    Clusters& clusters, Params& params, HyperParams& HP,
    InitMode init_mode, int pbwt_match_len,
    const SeqArray& x_train, int max_batch_size, int max_train_steps,
    Json& json
) {
    // Init clusters.
    auto t0 = std::chrono::steady_clock::now();
    switch (init_mode) {
        case InitMode::block:
            clusters.block_init(x_train);
            break;
        case InitMode::pbwt:
            if (HP.K != 2) { throw std::invalid_argument("PBWT init supports K=2 only."); }
            clusters.pbwt_init(x_train, pbwt_match_len);
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

    // Train.
    EarlyStopping early_stop{2, false, 1.0, max_train_steps};
    double elbo = 0.0;

    std::vector<Json> train_log;
    while (!early_stop.converged()) {
        auto t0 = std::chrono::steady_clock::now();
        expect_step(HP, params, clusters);
        auto t1 = std::chrono::steady_clock::now();
        max_step(x_train, clusters, params, HP, max_batch_size);
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
        double mean_nR = static_cast<double>(clusters.nR) / HP.L;
        std::cerr << early_stop.step << ": elbo=" << elbo << " mean_nR=" << mean_nR
            << " t_expect=" << t_expect << "ms t_max=" << t_max
            << "ms t_elbo=" << t_elbo << "ms t_step=" << t_step << "ms\n";
        train_log.emplace_back();
        train_log[train_log.size()-1].add("elbo", elbo).add("mean_nR", mean_nR)
            .add("t_max", t_max).add("t_expect", t_expect).add("t_elbo", t_elbo).add("t_step", t_step);
    }
    json.add("train_log", train_log);

    Json param_log;
    param_log.add("mu_alpha", params.mu_alpha).add("mu_gamma", params.mu_gamma).add("mu_d", params.mu_d)
        .add("alpha_eps", params.alpha_eps).add("beta_eps", params.beta_eps);
    json.add("params", param_log);
}


void get_viterbi_impute_probs(
    const SeqArray& x, int i, const std::unordered_map<int, int>& obs_ls,
    ViterbiBuffers& viterbi_bufs, std::vector<double>& seq_probs,
    const Clusters& clusters, const Params& params, const HyperParams& HP
) {
    get_viterbi_path(x, i, &obs_ls, viterbi_bufs, clusters, params, HP);
    int n_masked_ls = HP.L - obs_ls.size();
    int masked_l = 0;
    for (int l = 0; l < HP.L; ++l) {
        if (obs_ls.contains(l)) { continue; }
        for (int k = 0; k < HP.K; ++k) {
            double p = get_cluster_emission_ll(viterbi_bufs.path[2 * l], k, l, clusters, params, HP);
            seq_probs[idx2d(masked_l, k, HP.K)] = p;
        }
        ++masked_l;
    }
    normalize_ll(seq_probs, n_masked_ls, HP.K);
}

void impute(
    const SeqArray& x_val, const std::unordered_map<int, int>& obs_ls,
    bool viterbi,
    const Clusters& clusters, const Params& params, const HyperParams& HP,
    const std::string& prob_file
) {
    int n_masked_ls = HP.L - static_cast<int>(obs_ls.size());
    ImputeProbWriter prob_writer(prob_file, x_val.N, n_masked_ls);
    std::vector<double> seq_probs(n_masked_ls * HP.K);

    if (viterbi) {
        ViterbiBuffers viterbi_bufs(clusters.next_cluster_id, HP.L);
        for (int i = 0; i < x_val.N; ++i) {
            get_viterbi_impute_probs(x_val, i, obs_ls, viterbi_bufs, seq_probs, clusters, params, HP);
            prob_writer.write_row(seq_probs);
        }
    }
    else {
        FwdBkwdBuffers fwd_bkwd_bufs(clusters.next_cluster_id, HP.L);
        for (int i = 0; i < x_val.N; ++i) {
            fwd_bkwd(x_val, i, obs_ls, fwd_bkwd_bufs, seq_probs, clusters, params, HP);
            prob_writer.write_row(seq_probs);
        }
    }
    prob_writer.finish();
}


std::unordered_map<int, int> read_obs_ls(const std::string& fname, int n_obs_ls, int n_total_loci) {
    std::ifstream file(fname);
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
    // Parse args.
    CLI::App app{"DFCP imputation"};

    std::string ref_file;
    app.add_option("ref_file", ref_file, "Reference sequence file.")->required();

    std::optional<std::string> target_file;
    std::optional<std::string> obs_loci_file;
    std::optional<std::string> prob_file;
    app.add_option("target_file", target_file, "Target sequence file.");
    app.add_option("obs_loci_file", obs_loci_file, "Observed loci file.");
    app.add_option("prob_file", prob_file, "Imputation probability output.");

    HyperParams HP{.N=0, .L=0, .K=2};
    app.add_option("--tau_1", HP.tau_1, "Alpha prior shape.");
    app.add_option("--tau_2", HP.tau_2, "Alpha prior rate.");
    app.add_option("--v_1", HP.v_1, "Discount prior first shape.");
    app.add_option("--v_2", HP.v_2, "Discount prior second shape.");
    app.add_option("--phi_1", HP.phi_1, "Gamma prior shape.");
    app.add_option("--phi_2", HP.phi_2, "Gamma prior rate.");

    std::string mode = "hard";
    app.add_option("--mode", mode, "Emission mode: hard, noisy, or soft.")
        ->check(CLI::IsMember({"hard", "noisy", "soft"}));
    app.add_option("--lambda_1", HP.lambda_1, "Noisy emission prior first shape.");
    app.add_option("--lambda_2", HP.lambda_2, "Noisy emission prior second shape.");

    std::string init = "pbwt";
    int pbwt_match_len = 20;
    app.add_option("--init", init, "Initialization: viterbi, block, or pbwt.")
        ->check(CLI::IsMember({"viterbi", "block", "pbwt"}));
    app.add_option("--pbwt_match_len", pbwt_match_len, "PBWT match len.");

    int max_batch_size = 4;
    int max_train_steps = 3;
    app.add_option("--max_batch_size", max_batch_size, "Maximum sequence batch size.");
    app.add_option("--max_train_steps", max_train_steps, "Maximum training steps.");

    bool viterbi_impute = false;
    app.add_flag("--viterbi_impute", viterbi_impute, "Viterbi impute (default Forward-Backward).");

    std::optional<std::string> r_assign_file;
    app.add_option("--output_r_assign", r_assign_file, "Final R assignment file.");

    CLI11_PARSE(app, argc, argv);

    bool run_impute = target_file.has_value();
    if (run_impute != obs_loci_file.has_value() || run_impute != prob_file.has_value()) {
        throw std::invalid_argument("Imputation requires target, observed loci, and prob files.");
    }
    EmitMode emit_mode = mode == "hard" ? EmitMode::hard
        : mode == "noisy" ? EmitMode::noisy : EmitMode::soft;
    InitMode init_mode = init == "viterbi" ? InitMode::viterbi
        : init == "block" ? InitMode::block : InitMode::pbwt;

    // Read files.
    SeqArray x_train{read_seq_file(ref_file)};
    std::optional<SeqArray> x_val;
    std::unordered_map<int, int> obs_ls;
    if (run_impute) {
        x_val.emplace(read_seq_file(*target_file));
        obs_ls = read_obs_ls(*obs_loci_file, x_val->L, x_train.L);
    }

    Json json;
    json.add("ref_file", ref_file);

    std::cerr << "ref_file=" << ref_file << '\n';
    if (run_impute) {
        std::cerr << " target_file=" << *target_file
            << "\nobs_loci_file=" << *obs_loci_file
            << "\nprob_file=" << *prob_file << '\n';
        json.add("target_file", *target_file).add("obs_loci_file", *obs_loci_file);
    }
    HP.N = x_train.N;
    HP.L = x_train.L;
    std::cerr << HP << '\n';

    // Init params and clusters, train, impute.
    Params params{HP};
    Clusters clusters{HP, emit_mode};

    train_dfcp(
        clusters, params, HP,
        init_mode, pbwt_match_len,
        x_train, max_batch_size, max_train_steps,
        json
    );

    if (r_assign_file) {
        write_r_assign(*r_assign_file, clusters);
    }

    if (run_impute) {
        auto t0 = std::chrono::steady_clock::now();
        impute(*x_val, obs_ls, viterbi_impute, clusters, params, HP, *prob_file);
        auto t1 = std::chrono::steady_clock::now();
        auto t_impute = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        std::cerr << "t_impute=" << t_impute << "ms\n";
        json.add("t_impute", t_impute);
    }
    std::cout << json.str() << '\n';
}
