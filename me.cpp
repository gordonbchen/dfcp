#include <iostream>
#include <cassert>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <limits>
#include <iomanip>
#include <string_view>
#include <random>

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/logistic_sigmoid.hpp>
#include <boost/math/tools/minima.hpp>


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
};

std::ostream& operator<<(std::ostream& os, const HyperParams& HP) {
    return (
        os << "HP: N=" << HP.N << ", L=" << HP.L << ", K=" << HP.K
        << ", tau_1=" << HP.tau_1 << ", tau_2=" << HP.tau_2 << ", v_1=" << HP.v_1 << ", v_2=" << HP.v_2
        << ", phi_1=" << HP.phi_1 << ", phi_2=" << HP.phi_2
    );
}


struct Params {
    double mu_alpha;
    double sigma2_alpha;
    double mu_log_alpha;
    double sigma2_log_alpha;

    std::vector<double> mu_d;
    std::vector<double> sigma2_d;
    std::vector<double> mu_log_d;
    std::vector<double> sigma2_logit_d;

    std::vector<double> mu_gamma;
    std::vector<double> sigma2_gamma;
    std::vector<double> mu_log_gamma;
    std::vector<double> sigma2_log_gamma;

    Params(const HyperParams& HP) :
        mu_alpha(HP.tau_1 / HP.tau_2),
        sigma2_alpha(HP.tau_1 / (HP.tau_2*HP.tau_2)),
        mu_log_alpha(boost::math::digamma(HP.tau_1) - std::log(HP.tau_2)),
        sigma2_log_alpha(boost::math::trigamma(HP.tau_1)),

        mu_d(HP.L-1, HP.v_1 / (HP.v_1 + HP.v_2)),
        sigma2_d(HP.L-1, (HP.v_1*HP.v_2) / (std::pow(HP.v_1+HP.v_2, 2) * (HP.v_1+HP.v_2+1))),
        mu_log_d(HP.L-1, boost::math::digamma(HP.v_1) - boost::math::digamma(HP.v_1 + HP.v_2)),
        sigma2_logit_d(HP.L-1, boost::math::trigamma(HP.v_1) + boost::math::trigamma(HP.v_2)),

        mu_gamma(HP.L, HP.phi_1 / HP.phi_2),
        sigma2_gamma(HP.L, HP.phi_1 / std::pow(HP.phi_2, 2)),
        mu_log_gamma(HP.L, boost::math::digamma(HP.phi_1) - std::log(HP.phi_2)),
        sigma2_log_gamma(HP.L, boost::math::trigamma(HP.phi_1))
    {}
};


struct Cluster {
    std::unordered_set<int> seqs;
    const bool is_r;
    const int l;
    const int emission;

    std::unordered_set<Cluster*> parents;
    std::unordered_set<Cluster*> children;

    Cluster(std::unordered_set<int> seqs_, bool is_r_, int l_, int emission_) :
        seqs(std::move(seqs_)), is_r(is_r_), l(l_), emission(emission_)
    {
        assert(is_r == (emission != -1) && "only r cluster can have emissions.");
    }

    void add_child(Cluster *child) {
        children.insert(child);
        child->parents.insert(this);
    }
};


size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
}

double delta_Elogx(double mu, double sigma2, double a, double b) {
    double x = a*mu + b;
    return std::log(x) - 0.5*sigma2*a*a / (x*x);
}

struct Msg {
    double ll;
    Cluster* next;
};

double get_msg_ll(const std::unordered_map<Cluster*, Msg>& msgs, Cluster* c) {
    auto it = msgs.find(c);
    if (it == msgs.end()) {
        return -std::numeric_limits<double>::infinity();
    }
    return it->second.ll;
}

struct Clusters {
    std::unordered_map<Cluster*, std::unique_ptr<Cluster>> all_clusters;
    const HyperParams& HP;
    std::vector<char> modes;
    std::vector<int> nk;
    std::vector<Cluster*> r_assign;
    std::vector<Cluster*> q_assign;
    std::vector<std::unordered_set<Cluster*>> rs;
    std::vector<std::unordered_set<Cluster*>> qs;
    std::vector<std::unordered_set<Cluster*>> rs_by_emit;
    int nR;

    Clusters(const HyperParams& HP_, const std::vector<char>& x) :
        HP(HP_),
        modes(HP.L),
        nk(HP.L * HP.K, 0),
        r_assign(HP.N * HP.L, nullptr),
        q_assign(HP.N * (HP.L-1), nullptr),
        rs(HP.L), qs(HP.L-1),
        rs_by_emit(HP.L * HP.K),
        nR(0)
    {
        // Count modes.
        std::vector<int> counts(HP.K, 0);
        for (int l = 0; l < HP.L; ++l) {
            std::fill(counts.begin(), counts.end(), 0);
            for (int i = 0; i < HP.N; ++i) {
                int idx = idx2d(i, l, HP.L);
                if (x[idx] != -1) {
                    ++counts[x[idx]];
                }
            }
            auto max_it = std::max_element(counts.begin(), counts.end());
            assert(*max_it > 0 && "No valid alleles at loc.");
            modes[l] = std::distance(counts.begin(), max_it);
        }

        // Block init.
        std::unordered_set<int> seqs;
        seqs.reserve(HP.N);
        for (int i = 0; i < HP.N; ++i) {
            seqs.insert(i);
        }

        Cluster* r = create_cluster(seqs, true, 0, modes[0]);
        Cluster* q = nullptr;
        for (int l = 0; l < HP.L-1; ++l) {
            q = create_cluster(seqs, false, l, -1);
            r->add_child(q);

            r = create_cluster(seqs, true, l+1, modes[l+1]);
            q->add_child(r);
        }
    }

    Cluster* create_cluster(std::unordered_set<int> seqs, bool is_r, int l, int emission) {
        std::unique_ptr<Cluster> u_ptr = std::make_unique<Cluster>(std::move(seqs), is_r, l, emission);
        Cluster* ptr = u_ptr.get();
        all_clusters[ptr] = std::move(u_ptr);

        if (!is_r) {
            for (const int& i : ptr->seqs) {
                q_assign[idx2d(i, l, HP.L-1)] = ptr;
            }
            qs[l].insert(ptr);
            return ptr;
        }

        for (const int& i : ptr->seqs) {
            r_assign[idx2d(i, l, HP.L)] = ptr;
        }
        rs[l].insert(ptr);
        ++nk[idx2d(l, emission, HP.K)];
        rs_by_emit[idx2d(l, emission, HP.K)].insert(ptr);
        ++nR;
        return ptr;
    }

    void cluster_add(Cluster* cluster, int idx) {
        cluster->seqs.insert(idx);

        if (cluster->is_r) {
            r_assign[idx2d(idx, cluster->l, HP.L)] = cluster;
            return;
        }
        q_assign[idx2d(idx, cluster->l, HP.L-1)] = cluster;
    }

    void cluster_remove(Cluster* cluster, int idx) {
        cluster->seqs.erase(idx);

        if (cluster->is_r) {
            r_assign[idx2d(idx, cluster->l, HP.L)] = nullptr;
        }
        else {
            q_assign[idx2d(idx, cluster->l, HP.L-1)] = nullptr;
        }

        if (cluster->seqs.size() > 0) {
            return;
        }

        // Delete cluster.
        for (Cluster* parent : cluster->parents) {
            parent->children.erase(cluster);
        }
        for (Cluster* child: cluster->children) {
            child->parents.erase(cluster);
        }

        if (cluster->is_r) {
            rs[cluster->l].erase(cluster);
            --nk[idx2d(cluster->l, cluster->emission, HP.K)];
            rs_by_emit[idx2d(cluster->l, cluster->emission, HP.K)].erase(cluster);
            --nR;
        }
        else {
            qs[cluster->l].erase(cluster);
        }

        all_clusters.erase(cluster);
    }

    void max_step(const std::vector<char>& x, const Params& params) {
        for (int i = 0; i < HP.N; ++i) {
            for (int l = 0; l < HP.L; ++l) {
                cluster_remove(r_assign[idx2d(i, l, HP.L)], i);
                if (l == HP.L-1) {
                    break;
                }
                cluster_remove(q_assign[idx2d(i, l, HP.L-1)], i);
            }
            viterbi_seq(x.begin() + idx2d(i, 0, HP.L), i, params);
        }
    }

    void viterbi_seq(std::vector<char>::const_iterator xi, int i, const Params& params) {
        std::vector<std::unordered_map<Cluster*, Msg>> a_msgs(HP.L);
        std::vector<std::unordered_map<Cluster*, Msg>> b_msgs(HP.L-1);
        for (int l = HP.L-1; l >= 0; --l) {
            // Likelihood for new cluster.
            auto& ma = a_msgs[l];
            int emission = xi[l] == -1 ? modes[l] : xi[l];
            double new_a_ll = (
                delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, nk[idx2d(l, emission, HP.K)])
                - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, rs[l].size())
            );

            std::unordered_set<Cluster*>& matching_as = (
                xi[l] == -1 ? rs[l] : rs_by_emit[idx2d(l, emission, HP.K)]
            );
            if (l == HP.L-1) {
                ma[nullptr] = Msg{new_a_ll, nullptr};
                for (Cluster *a : matching_as) {
                    ma[a] = Msg{0.0, nullptr};
                }
                continue;
            }

            // b messages.
            auto& next_ma = a_msgs[l+1];
            auto& mb = b_msgs[l];
            for (Cluster* b : qs[l]) {
                assert(b->children.size() == 1 && "b clusters should only have 1 child.");
                Cluster* next_a = *b->children.begin();
                mb[b] = Msg{get_msg_ll(next_ma, next_a), next_a};
             }

            int nQl = qs[l].size();
            double mu_y = params.mu_alpha + nQl*params.mu_d[l];
            double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
            double elogy = delta_Elogx(mu_y, sigma2_y, 1.0, 0.0);

            Cluster* best_a = nullptr;
            double best_a_ll = params.mu_log_alpha + next_ma.at(nullptr).ll;
            std::unordered_set<Cluster*>& matching_next_as = (
                xi[l+1] == -1 ? rs[l+1] : rs_by_emit[idx2d(l+1, xi[l+1], HP.K)]
            );
            for (Cluster *a : matching_next_as) {
                double nCl = a->parents.size();
                double ll = params.mu_log_d[l] + std::log(nCl) + get_msg_ll(next_ma, a);
                if (ll > best_a_ll) {
                    best_a = a;
                    best_a_ll = ll;
                }
            }
            double new_b_ll = -elogy + best_a_ll;
            mb[nullptr] = Msg{new_b_ll, best_a};

            // a messages.
            ma[nullptr] = Msg{new_a_ll + new_b_ll, nullptr};
            for (Cluster* a : matching_as) {
                Cluster* best_b = nullptr;
                double nFl = a->children.size();
                double best_b_ll = std::log(nFl) + params.mu_log_d[l] + mb[nullptr].ll;
                for (Cluster* b : a->children) {
                    double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->seqs.size()) + get_msg_ll(mb, b);
                    if (ll > best_b_ll) {
                        best_b = b;
                        best_b_ll = ll;
                    }
                }
                best_b_ll -= std::log(static_cast<double>(a->seqs.size()));
                ma[a] = Msg{best_b_ll, best_b};
            }
        }

        // Viterbi path.
        Cluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
            [](const auto& a, const auto& b) { return a.second.ll < b.second.ll; }
        )->first;
        Cluster* b = nullptr;
        Cluster* next_a = nullptr;

        for (int l = 0; l < HP.L-1; ++l) {
            b = a_msgs[l].at(a).next;
            if (a == nullptr) {
                a = next_a;
            }
            next_a = b_msgs[l].at(b).next;

            if (b == nullptr) {
                b = create_cluster(std::unordered_set<int>{i}, false, l, -1);
            }
            else {
                cluster_add(b, i);
            }

            if (l == 0 && a == nullptr) {
                int emission = xi[l] == -1 ? modes[l] : xi[l];
                a = create_cluster(std::unordered_set<int>{i}, true, l, emission);
            }
            else {
                cluster_add(a, i);
            }
            a->add_child(b);

            a = next_a;
            if (next_a == nullptr) {
                int next_emission = xi[l+1] == -1 ? modes[l+1] : xi[l+1];
                next_a = create_cluster(std::unordered_set<int>{i}, true, l+1, next_emission);
            }
            else {
                cluster_add(next_a, i);
            }
            b->add_child(next_a);
        }
    }

    void add_seqs(const std::vector<char>& x, const Params& params, HyperParams& HP) {
        int n = x.size() / HP.L;
        int old_N = HP.N;
        HP.N += n;
        r_assign.resize(r_assign.size() + x.size(), nullptr);
        q_assign.resize(q_assign.size() + n * (HP.L-1), nullptr);

        for (int i = 0; i < n; ++i) {
            viterbi_seq(x.begin() + i*HP.L, old_N + i, params);
        }
    }
};


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


template <typename F_nll, typename F_d2>
void laplace_log_approx(
    F_nll nll_log_func, F_d2 ll_log_d2_func,
    double& mu, double& sigma2,
    double& logx_mode, double& logx_var
) {
    // TODO: bits and max iter for find min.
    logx_mode = boost::math::tools::brent_find_minima(nll_log_func, -10.0, 10.0, 30).first;
    logx_var = -1.0 / ll_log_d2_func(std::exp(logx_mode));
    assert(logx_var > 0.0);

    mu = std::exp(logx_mode + 0.5*logx_var);
    sigma2 = (std::exp(logx_var) - 1.0) * std::exp(2.0*logx_mode + logx_var);
    assert(sigma2 > 0.0);
}

double delta_ElogGamma_invx(double mu, double sigma2, double a, double b) {
    double x = a/mu + b;
    double d2 = boost::math::trigamma(x) * (a*a) / std::pow(mu, 4);
    d2 += boost::math::digamma(x) * 2*a / std::pow(mu, 3);
    return std::lgamma(x) + 0.5*sigma2 * d2;
}

double delta_ElogGamma_invx_d2_x(double mu, double sigma2, double a, double b) {
    double x = a/mu + b;
    double d2 = boost::math::trigamma(x) / (mu*mu);
    // TODO: polygamma vs splitting into sum and delta approxs.
    d2 += 0.5*sigma2 * (
        boost::math::polygamma(3, x) * (a*a) / std::pow(mu, 6)
        + boost::math::polygamma(2, x) * (6*a) / std::pow(mu, 5)
        + boost::math::polygamma(1, x) * 6 / std::pow(mu, 4)
    );
    return d2;
}

double ll_log_alpha(double log_alpha, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    double alpha = std::exp(log_alpha);

    double ll = std::lgamma(alpha) - std::lgamma(alpha + HP.N);
    ll += (clusters.nR + HP.tau_1 - 1)*std::log(alpha) - HP.tau_2*alpha;
    for (int l = 0; l < HP.L-1; ++l) {
        ll += delta_ElogGamma_invx(params.mu_d[l], params.sigma2_d[l], alpha, 0.0);
        ll -= delta_ElogGamma_invx(params.mu_d[l], params.sigma2_d[l], alpha, clusters.qs[l].size());
    }

    return ll + log_alpha;
}

double ll_log_alpha_d2(double alpha, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    double d2 = boost::math::trigamma(alpha) - boost::math::trigamma(alpha + HP.N);
    d2 += (1 - HP.tau_1 - clusters.nR) / (alpha*alpha);
    for (int l = 0; l < HP.L-1; ++l) {
        d2 += delta_ElogGamma_invx_d2_x(params.mu_d[l], params.sigma2_d[l], alpha, 0.0);
        d2 -= delta_ElogGamma_invx_d2_x(params.mu_d[l], params.sigma2_d[l], alpha, clusters.qs[l].size());
    }
    return alpha*alpha * d2 - 1;
}

double ll_log_gammal(double log_gamma, int l, const HyperParams& HP, const Clusters& clusters) {
    double gamma = std::exp(log_gamma);

    double ll = (HP.phi_1-1)*std::log(gamma) - HP.phi_2*gamma;
    ll += std::lgamma(HP.K*gamma) - std::lgamma(HP.K*gamma + clusters.rs[l].size());
    for (int k = 0; k < HP.K; ++k) {
        ll += std::lgamma(gamma + clusters.nk[idx2d(l, k, HP.K)]);
    }
    ll -= HP.K * std::lgamma(gamma);

    return ll + log_gamma;
}

double ll_log_gammal_d2(double gamma, int l, const HyperParams& HP, const Clusters& clusters) {
    double d2 = (1.0 - HP.phi_1) / (gamma*gamma);
    d2 += HP.K*HP.K * (boost::math::trigamma(HP.K*gamma) - boost::math::trigamma(HP.K*gamma + clusters.rs[l].size()));
    for (int k = 0; k < HP.K; ++k) {
        d2 += boost::math::trigamma(gamma + clusters.nk[idx2d(l, k, HP.K)]);
        d2 -= boost::math::trigamma(gamma);
    }
    return gamma*gamma * d2 - 1;
}

double delta_ElogGamma_x(double mu, double sigma2, double a = 1.0, double b = 0.0) {
    double x = a*mu + b;
    return std::lgamma(x) + 0.5 * sigma2 * a * a * boost::math::trigamma(x);
}

double delta_ElogGamma_x_d2_invx(double mu, double sigma2, double a, double b) {
    double x = a*mu + b;
    double d2 = boost::math::trigamma(x) * mu*mu * std::pow(a, 4);
    d2 += boost::math::digamma(x) * 2*mu * std::pow(a, 3);
    d2 += 0.5 * sigma2 * (
        boost::math::polygamma(3, x) * mu*mu * std::pow(a, 6)
        + boost::math::polygamma(2, x) * 6*mu * std::pow(a, 5)
        + boost::math::trigamma(x) * 6 * std::pow(a, 4)
    );
    return d2;
}

double ll_logit_dl(double logit_dl, int l, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    double d = boost::math::logistic_sigmoid(logit_dl);

    int nQl = clusters.qs[l].size();
    int nRl1 = clusters.rs[l].size() + clusters.rs[l+1].size();
    double ll = (nQl - nRl1 + HP.v_1 - 1.0) * std::log(d);
    ll += (HP.v_2 - 1.0) * std::log(1.0 - d);
    ll -= nQl * std::lgamma(1.0 - d);
    for (Cluster* b : clusters.qs[l]) {
        ll += std::lgamma(b->seqs.size() - d);
    }
    ll += delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0/d, 0.0);
    ll -= delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0/d, nQl);

    return ll + std::log(d) + std::log(1.0-d);
}

double ll_logit_dl_d2(double d, int l, const HyperParams& HP, const Params& params, const Clusters& clusters) {
    int nQl = clusters.qs[l].size();
    int nRl1 = clusters.rs[l].size() + clusters.rs[l+1].size();
    double d2 = (nRl1 - nQl + 1.0 - HP.v_1) / (d*d);
    d2 -= (HP.v_2 - 1.0) / std::pow(1.0 - d, 2);
    d2 -= nQl * boost::math::trigamma(1.0 - d);
    for (Cluster* b : clusters.qs[l]) {
        d2 += boost::math::trigamma(b->seqs.size() - d);
    }

    d2 += delta_ElogGamma_x_d2_invx(params.mu_alpha, params.sigma2_alpha, 1.0/d, 0.0);
    d2 -= delta_ElogGamma_x_d2_invx(params.mu_alpha, params.sigma2_alpha, 1.0/d, nQl);

    d2 *= std::pow(d * (1.0-d), 2);
    d2 -= std::pow(1.0 - 2.0*d, 2) + 2.0*d*(1.0-d);
    return d2;
}

void expect_step(const HyperParams& HP, Params& params, const Clusters& clusters) {
    // alpha update.
    laplace_log_approx(
        [&](double log_alpha) -> double { return -ll_log_alpha(log_alpha, HP, params, clusters); },
        [&](double alpha) -> double { return ll_log_alpha_d2(alpha, HP, params, clusters); },
        params.mu_alpha, params.sigma2_alpha,
        params.mu_log_alpha, params.sigma2_log_alpha
    );

    for (int l = 0; l < HP.L; ++l) {
        // gamma_l update.
        laplace_log_approx(
            [&](double log_gammal) { return -ll_log_gammal(log_gammal, l, HP, clusters); },
            [&](double gammal) { return ll_log_gammal_d2(gammal, l, HP, clusters); },
            params.mu_gamma[l], params.sigma2_gamma[l],
            params.mu_log_gamma[l], params.sigma2_log_gamma[l]
        );
        if (l >= HP.L-1) {
            break;
        }

        // d_l update in logit space.
        double logit_dl_mode = boost::math::tools::brent_find_minima(
            [&](double logit_dl) {return -ll_logit_dl(logit_dl, l, HP, params, clusters); },
            -10.0, 10.0, 30
        ).first;
        double dl_mode = boost::math::logistic_sigmoid(logit_dl_mode);
        params.sigma2_logit_d[l] = -1.0 / ll_logit_dl_d2(dl_mode, l, HP, params, clusters);
        assert(params.sigma2_logit_d[l] > 0.0);

        params.mu_d[l] = dl_mode + 0.5*params.sigma2_logit_d[l] * (1.0 - 2.0*dl_mode)*dl_mode*(1.0 - dl_mode);
        params.sigma2_d[l] = (
            dl_mode*dl_mode
            + 0.5*params.sigma2_logit_d[l] * (4.0*dl_mode - 6.0*dl_mode*dl_mode) * dl_mode * (1.0 - dl_mode)
            - params.mu_d[l]*params.mu_d[l]
        );
        assert(params.sigma2_d[l] > 0.0);
        params.mu_log_d[l] = delta_Elogx(params.mu_d[l], params.sigma2_d[l], 1.0, 0.0);
    }
}


double normal_entropy(double sigma2) {
    assert(sigma2 > 0.0);
    constexpr double pi = 3.141592653589793238462643383279502884;
    return 0.5 * std::log(2.0 * pi * std::exp(1.0) * sigma2);
}

double betaln(double a, double b) {
    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

double calc_elbo(const HyperParams& HP, const Params& params, const Clusters& clusters) {
    // alpha.
    double elbo = delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha);
    elbo -= delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0, HP.N);
    elbo += (clusters.nR + HP.tau_1) * params.mu_log_alpha;
    elbo -= HP.tau_2 * params.mu_alpha;
    elbo += HP.tau_1 * std::log(HP.tau_2) - std::lgamma(HP.tau_1);

    for (int l = 0; l < HP.L - 1; ++l) {
        // d.
        // -1s cancel out with d_l entropy.
        int nQl = clusters.qs[l].size();
        int nRl1 = clusters.rs[l].size() + clusters.rs[l+1].size();
        elbo += (nQl - nRl1 + HP.v_1) * params.mu_log_d[l];
        elbo += HP.v_2 * delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1.0, 1.0);
        elbo -= betaln(HP.v_1, HP.v_2);
        elbo -= nQl * delta_ElogGamma_x(params.mu_d[l], params.sigma2_d[l], -1.0, 1.0);

        for (Cluster* b : clusters.qs[l]) {
            elbo += delta_ElogGamma_x(params.mu_d[l], params.sigma2_d[l], -1.0, b->seqs.size());
        }

        // alpha and d term.
        double z = params.mu_alpha / params.mu_d[l];
        double dd2 = boost::math::digamma(z) * (2.0 * z / (params.mu_d[l] * params.mu_d[l]));
        dd2 += boost::math::trigamma(z) * (z * z / (params.mu_d[l] * params.mu_d[l]));
        elbo += delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0 / params.mu_d[l], 0.0);
        elbo += 0.5 * params.sigma2_d[l] * dd2;

        z += nQl;
        dd2 = boost::math::digamma(z) * (2.0 * params.mu_alpha / std::pow(params.mu_d[l], 3));
        dd2 += boost::math::trigamma(z) * (params.mu_alpha * params.mu_alpha / std::pow(params.mu_d[l], 4));
        elbo -= delta_ElogGamma_x(params.mu_alpha, params.sigma2_alpha, 1.0 / params.mu_d[l], nQl);
        elbo -= 0.5 * params.sigma2_d[l] * dd2;
    }

    // gamma.
    elbo += HP.L * (HP.phi_1 * std::log(HP.phi_2) - std::lgamma(HP.phi_1));
    for (int l = 0; l < HP.L; ++l) {
        elbo += HP.phi_1 * params.mu_log_gamma[l];
        elbo -= HP.phi_2 * params.mu_gamma[l];
        elbo += delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, 0.0);
        elbo -= delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, clusters.rs[l].size());
        for (int k = 0; k < HP.K; ++k) {
            elbo += delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, clusters.nk[idx2d(l, k, HP.K)]);
        }
        elbo -= HP.K * delta_ElogGamma_x(params.mu_gamma[l], params.sigma2_gamma[l]);
    }

    // Clusters.
    for (Cluster* a : clusters.rs[0]) {
        elbo += std::lgamma(static_cast<double>(a->seqs.size()));
    }
    for (int l = 0; l < HP.L - 1; ++l) {
        for (Cluster* a : clusters.rs[l]) {
            elbo += std::lgamma(static_cast<double>(a->children.size()));
            elbo -= std::lgamma(static_cast<double>(a->seqs.size()));
        }
        for (Cluster* a : clusters.rs[l + 1]) {
            elbo += std::lgamma(static_cast<double>(a->parents.size()));
        }
    }

    // Variational entropy.
    elbo += normal_entropy(params.sigma2_log_alpha);
    for (int l = 0; l < HP.L; ++l) {
        elbo += normal_entropy(params.sigma2_log_gamma[l]);
    }
    for (int l = 0; l < HP.L - 1; ++l) {
        elbo += normal_entropy(params.sigma2_logit_d[l]);
    }
    return elbo;
}


struct SparseX {
    size_t i;
    size_t l;
    char x;
    SparseX(size_t i_, size_t l_, char x_) : i(i_), l(l_), x(x_) {}
};

int main(int argc, char *argv[]) {
    // Read seq file.
    assert(argc >= 2 && "Requires sequence file.");
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "seq_file=" << argv[1] << '\n';
    std::ifstream seq_file(argv[1]);
    assert(seq_file.is_open() && "Failed to open sequence file.");

    int N = 0;
    int L = 0;
    std::vector<char> x;
    std::string line;
    while (std::getline(seq_file, line)) {
        if (N == 0) {
            L = line.length();
        }
        for (char c : line) {
            x.push_back(c - '0');
        }
        ++N;
    }
    int K = *std::max_element(x.begin(), x.end()) + 1;

    // Read hyperparams.
    HyperParams HP{.N=N, .L=L, .K=K};
    // TODO: mask and imputation check.
    double val = -1.0;
    double mask = -1.0;
    std::unordered_map<std::string_view, double*> args = {
        {"--tau_1", &HP.tau_1}, {"--tau_2", &HP.tau_2},
        {"--v_1", &HP.v_1}, {"--v_2", &HP.v_2},
        {"--phi_1", &HP.phi_1}, {"--phi_2", &HP.phi_2},
        {"--val", &val}, {"--mask", &mask}
    };

    char* end_ptr = nullptr;
    int i = 2;
    while (i < argc) {
        auto it = args.find(argv[i]);
        assert(it != args.end() && "Arg not recognized.");

        assert((i+1 < argc) && "Arg has no value.");
        *it->second = std::strtod(argv[i+1], &end_ptr);
        assert((end_ptr != argv[i+1]) && "Failed to parse arg value double.");

        args.erase(it);
        i += 2;
    }

    std::cout << HP << '\n';
    std::cout << "val=" << val << ", mask=" << mask << '\n';

    // Split val for imputation.
    bool do_val = val > 0.0;
    assert((do_val == (mask > 0.0)) && "If imputation validation frac is positive, mask frac must be positive.");
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
        clusters.max_step(x_train, params);
        expect_step(HP, params, clusters);
        elbo = calc_elbo(HP, params, clusters);
        early_stop.update(elbo);
        std::cout << early_stop.step << ": " << "elbo=" << elbo << '\n';
    }

    // Impute.
    if (n_val_seqs > 0) {
        clusters.add_seqs(x_val_masked, params, HP);
        int n_correct = 0;
        for (SparseX& s : x_val_true) {
            if (s.x == clusters.r_assign[idx2d(n_train + s.i, s.l, HP.L)]->emission) {
                ++n_correct;
            }
        }
        std::cout << "Imputation acc: " << n_correct << '/' << n_masked_alleles << " = "
            << static_cast<double>(n_correct) / static_cast<double>(n_masked_alleles) << '\n';
    }
    return 0;
}

