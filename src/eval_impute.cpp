int get_max_prob_emission(std::vector<double>::const_iterator probs, int K) {
    int best_k = 0;
    double best_prob = probs[0];
    for (int k = 1; k < K; ++k) {
        if (probs[k] > best_prob) {
            best_k = k;
            best_prob = probs[k];
        }
    }
    return best_k;
}

struct R2Calc {
    int n = 0;
    double mean_x = 0.0;
    double mean_y = 0.0;
    double var_x = 0.0;
    double var_y = 0.0;
    double cov = 0.0;

    void update(double x, double y) {
        ++n;
        double delta_x = x - mean_x;
        double delta_y = y - mean_y;
        mean_x += delta_x / n;
        mean_y += delta_y / n;
        var_x += delta_x * (x - mean_x);
        var_y += delta_y * (y - mean_y);
        cov += delta_x * (y - mean_y);
    }

    double r2() const {
        if ((var_x <= 0.0) || (var_y <= 0.0)) {
            return -1.0;
        }
        return std::clamp((cov * cov) / (var_x * var_y), 0.0, 1.0);
    }
};

void impute(
    const Clusters& clusters, const Params& params, const HyperParams& HP,

    const std::vector<int8_t>& x_train,
    const std::vector<int8_t>& x_val, int n_val_seqs,
    const std::vector<int8_t>& x_val_true, const std::vector<int>& masked_ls,

    Json& json
) {
    std::vector<size_t> emission_counts{count_emissions(x_train, HP.N, HP.L, HP.K)};
    std::vector<int8_t> modes{get_emission_modes(emission_counts, HP.L, HP.K)};

    int n_viterbi_correct = 0;
    int n_fwd_bkwd_correct = 0;
    int n_mode_correct = 0;

    int n_masked_ls = masked_ls.size();
    int n_r2 = (HP.K == 2) ? n_masked_ls : 0;
    std::vector<R2Calc> fwd_bkwd_correlations(n_r2);
    std::vector<R2Calc> viterbi_correlations(n_r2);
    auto [minor_alleles, minor_allele_counts] = count_minor_alleles(emission_counts, masked_ls, n_masked_ls, HP.K);

    std::chrono::steady_clock::duration viterbi_impute_dur{};
    std::chrono::steady_clock::duration fwd_bkwd_impute_dur{};
    for (int i = 0; i < n_val_seqs; ++i) {
        auto t0 = std::chrono::steady_clock::now();
        std::vector<double> fwd_bkwd_probs{fwd_bkwd(x_val.begin() + i*HP.L, masked_ls, clusters, params, HP)};
        fwd_bkwd_impute_dur += std::chrono::steady_clock::now() - t0;

        t0 = std::chrono::steady_clock::now();
        std::vector<Cluster*> viterbi_clusters{get_viterbi_clusters(clusters, x_val.begin() + i*HP.L, HP, params)};
        std::vector<double> viterbi_probs;
        viterbi_probs.reserve(n_masked_ls * HP.K);
        for (int il = 0; il < n_masked_ls; ++il) {
            int l = masked_ls[il];
            for (int k = 0; k < HP.K; ++k) {
                double p = get_cluster_emission_ll(viterbi_clusters[2*l], k, l, clusters, params, HP);
                viterbi_probs.emplace_back(p);
            }
        }
        normalize_ll(viterbi_probs, n_masked_ls, HP.K);
        viterbi_impute_dur += std::chrono::steady_clock::now() - t0;

        for (int il = 0; il < n_masked_ls; ++il) {
            int l = masked_ls[il];
            int x_true = x_val_true[idx2d(i,il,n_masked_ls)];

            // Accuracy.
            if (x_true == get_max_prob_emission(viterbi_probs.begin() + il*HP.K, HP.K)) {
                ++n_viterbi_correct;
            }
            if (x_true == get_max_prob_emission(fwd_bkwd_probs.begin() + il*HP.K, HP.K)) {
                ++n_fwd_bkwd_correct;
            }
            if (x_true == modes[l]) {
                ++n_mode_correct;
            }

            // r^2 b/t minor allele and p(minor allele).
            if (HP.K == 2) {
                int minor_allele = minor_alleles[il];
                int is_minor_allele = x_true == minor_allele;
                fwd_bkwd_correlations[il].update(is_minor_allele, fwd_bkwd_probs[idx2d(il,minor_allele, HP.K)]);
                viterbi_correlations[il].update(is_minor_allele, viterbi_probs[idx2d(il,minor_allele, HP.K)]);
            }
        }
    }
    double n_masked_alleles = n_val_seqs * n_masked_ls;
    double viterbi_impute_acc = n_viterbi_correct / n_masked_alleles;
    double mode_impute_acc = n_mode_correct / n_masked_alleles;
    double fwd_bkwd_impute_acc = n_fwd_bkwd_correct / n_masked_alleles;

    std::vector<double> fwd_bkwd_r2s(n_r2, -1.0);
    std::vector<double> viterbi_r2s(n_r2, -1.0);
    double fwd_bkwd_r2_mean = 0.0;
    double viterbi_r2_mean = 0.0;
    int n_fwd_bkwd_r2_locs = 0;
    int n_viterbi_r2_locs = 0;
    if (HP.K == 2) {
        for (int il = 0; il < n_masked_ls; ++il) {
            fwd_bkwd_r2s[il] = fwd_bkwd_correlations[il].r2();
            viterbi_r2s[il] = viterbi_correlations[il].r2();
            if (fwd_bkwd_r2s[il] >= 0.0) {
                fwd_bkwd_r2_mean += fwd_bkwd_r2s[il];
                ++n_fwd_bkwd_r2_locs;
            }
            if (viterbi_r2s[il] >= 0.0) {
                viterbi_r2_mean += viterbi_r2s[il];
                ++n_viterbi_r2_locs;
            }
        }
        fwd_bkwd_r2_mean = (n_fwd_bkwd_r2_locs == 0) ? -1.0 : fwd_bkwd_r2_mean / n_fwd_bkwd_r2_locs;
        viterbi_r2_mean = (n_viterbi_r2_locs == 0) ? -1.0 : viterbi_r2_mean / n_viterbi_r2_locs;
        std::cerr << "fwd_bkwd_r2_mean=" << fwd_bkwd_r2_mean << " viterbi_r2_mean=" << viterbi_r2_mean << '\n';
        json.add("fwd_bkwd_r2s", fwd_bkwd_r2s).add("viterbi_r2s", viterbi_r2s)
            .add("fwd_bkwd_r2_mean", fwd_bkwd_r2_mean).add("viterbi_r2_mean", viterbi_r2_mean)
            .add("minor_allele_counts", minor_allele_counts);
    }

    auto t_viterbi_impute = std::chrono::duration_cast<std::chrono::milliseconds>(viterbi_impute_dur).count();
    auto t_fwd_bkwd_impute = std::chrono::duration_cast<std::chrono::milliseconds>(fwd_bkwd_impute_dur).count();
    std::cerr << "viterbi_impute_acc=" << viterbi_impute_acc << " t_viterbi_impute=" << t_viterbi_impute << "ms\n"
        << "fwd_bkwd_impute_acc=" << fwd_bkwd_impute_acc << " t_fwd_bkwd_impute=" << t_fwd_bkwd_impute << "ms\n"
        << "mode_impute_acc=" << mode_impute_acc << '\n';
    json.add("viterbi_impute_acc", viterbi_impute_acc).add("t_viterbi_impute", t_viterbi_impute)
        .add("fwd_bkwd_impute_acc", fwd_bkwd_impute_acc).add("t_fwd_bkwd_impute", t_fwd_bkwd_impute)
        .add("mode_impute_acc", mode_impute_acc);
}

