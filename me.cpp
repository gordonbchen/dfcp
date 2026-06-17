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

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>


struct HyperParams {
    int N;
    int L;
    int K;
    double tau_1;
    double tau_2;
    double v_1;
    double v_2;
    double phi_1;
    double phi_2;
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
    bool is_r;
    int l;
    int emission;

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


inline size_t idx2d(size_t r, size_t c, size_t width) {
    return r*width + c;
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

    Clusters(const HyperParams& HP_, const std::vector<char>& x) :
        HP(HP_),
        modes(HP.L),
        nk(HP.L * HP.K, 0),
        r_assign(HP.N * HP.L, nullptr),
        q_assign(HP.N * (HP.L-1), nullptr),
        rs(HP.L), qs(HP.L-1),
        rs_by_emit(HP.L * HP.K)
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
        }
        else {
            qs[cluster->l].erase(cluster);
        }

        all_clusters.erase(cluster);
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

        bool converged() {
            return steps_since_min > patience;
        }
};


inline double delta_Elogx(double mu, double sigma2, double a = 1.0, double b = 0.0) {
    double x = a*mu + b;
    return std::log(x) - 0.5*sigma2*a*a / (x*x);
}

struct Msg {
    double ll;
    Cluster* next;
};

double get_msg_ll(std::unordered_map<Cluster*, Msg>& msgs, Cluster* c) {
    auto it = msgs.find(c);
    if (it == msgs.end()) {
        return -std::numeric_limits<double>::infinity();
    }
    return it->second.ll;
}

void max_step(std::vector<char>& x, HyperParams& HP, Params& params, Clusters& clusters) {
    for (int i = 0; i < HP.N; ++i) {
        for (int l = 0; l < HP.L; ++l) {
            clusters.cluster_remove(clusters.r_assign[idx2d(i, l, HP.L)], i);
            if (l == HP.L-1) {
                break;
            }
            clusters.cluster_remove(clusters.q_assign[idx2d(i, l, HP.L-1)], i);
        }

        std::vector<std::unordered_map<Cluster*, Msg>> a_msgs(HP.L);
        std::vector<std::unordered_map<Cluster*, Msg>> b_msgs(HP.L-1);
        for (int l = HP.L-1; l >= 0; --l) {
            // Likelihood for new cluster.
            auto& ma = a_msgs[l];
            int emission = x[idx2d(i, l, HP.L)] == -1 ? clusters.modes[l] : x[idx2d(i, l, HP.L)];
            double new_a_ll = (
                delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], 1.0, clusters.nk[idx2d(l, emission, HP.K)])
                - delta_Elogx(params.mu_gamma[l], params.sigma2_gamma[l], HP.K, clusters.rs[l].size())
            );

            std::unordered_set<Cluster*>& matching_as = (
                x[idx2d(i, l, HP.L)] == -1 ? clusters.rs[l] : clusters.rs_by_emit[idx2d(l, emission, HP.K)]
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
            for (Cluster* b : clusters.qs[l]) {
                assert(b->children.size() == 1 && "b clusters should only have 1 child.");
                Cluster* next_a = *b->children.begin();
                mb[b] = Msg{get_msg_ll(next_ma, next_a), next_a};
             }

            int nQl = clusters.qs[l].size();
            double mu_y = params.mu_alpha + nQl*params.mu_d[l];
            double sigma2_y = params.sigma2_alpha + nQl*nQl * params.sigma2_d[l];
            double elogy = delta_Elogx(mu_y, sigma2_y);

            Cluster* best_a = nullptr;
            double best_a_ll = params.mu_log_alpha + next_ma.at(nullptr).ll;
            std::unordered_set<Cluster*>& matching_next_as = (
                x[idx2d(i, l+1, HP.L)] == -1 ?
                clusters.rs[l+1] : clusters.rs_by_emit[idx2d(l+1, x[idx2d(i, l+1, HP.L)], HP.K)]
            );
            for (Cluster *a : matching_next_as) {
                double ll = params.mu_log_d[l] + std::log(a->parents.size()) + get_msg_ll(next_ma, a);
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
                double best_b_ll = std::log(a->children.size()) + params.mu_log_d[l] + mb[nullptr].ll;
                for (Cluster* b : a->children) {
                    double ll = delta_Elogx(params.mu_d[l], params.sigma2_d[l], -1, b->seqs.size()) + get_msg_ll(mb, b);
                    if (ll > best_b_ll) {
                        best_b = b;
                        best_b_ll = ll;
                    }
                }
                ma[a] = Msg{-std::log(a->seqs.size()) + best_b_ll, best_b};
            }
        }

        // Viterbi path.
        Cluster* a = std::max_element(a_msgs[0].begin(), a_msgs[0].end(),
            [](const auto& a, const auto& b) {
                return a.second.ll < b.second.ll;
            }
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
                b = clusters.create_cluster(std::unordered_set<int>{i}, false, l, -1);
            }
            else {
                clusters.cluster_add(b, i);
            }

            if (l == 0 && a == nullptr) {
                int emission = x[idx2d(i, l, HP.L)] == -1 ? clusters.modes[l] : x[idx2d(i, l, HP.L)];
                a = clusters.create_cluster(std::unordered_set<int>{i}, true, l, emission);
            }
            else {
                clusters.cluster_add(a, i);
            }
            a->add_child(b);

            a = next_a;
            if (next_a == nullptr) {
                int next_emission = x[idx2d(i, l+1, HP.L)] == -1 ? clusters.modes[l+1] : x[idx2d(i, l+1, HP.L)];
                next_a = clusters.create_cluster(std::unordered_set<int>{i}, true, l+1, next_emission);
            }
            else {
                clusters.cluster_add(next_a, i);
            }
            b->add_child(next_a);
        }
    }
}


void expect_step(HyperParams& HP, Params& params, Clusters& clusters) {
}


double calc_elbo(HyperParams& HP, Params& params, Clusters& clusters) {
    return 0.0;
}


int main(int argc, char *argv[]) {
    // Read seq file.
    assert(argc == 2 && "Requires sequence file.");
    std::cout << "seq_file: " << argv[1] << '\n';
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
    HyperParams HP{.N=N, .L=L, .K=K, .tau_1=1.0, .tau_2=1.0, .v_1=1.0, .v_2=1.0, .phi_1=2.0, .phi_2=2.0};
    std::cout << HP << '\n';

    // Init params and clusters.
    Params params{HP};
    Clusters clusters{HP, x};

    EarlyStopping early_stop{3, false, 1e-5};
    double elbo = 0.0;
    while (!early_stop.converged()) {
        max_step(x, HP, params, clusters);
        expect_step(HP, params, clusters);
        elbo = calc_elbo(HP, params, clusters);
        early_stop.update(elbo);
        std::cout << early_stop.step << ": " << "elbo=" << elbo << '\n';
    }

    return 0;
}

