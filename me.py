from __future__ import annotations
import uuid
from argparse import ArgumentParser

import graphviz
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, special, optimize
from line_profiler import profile

from generate import generate, HyperParams


class Cluster:
    def __init__(
        self, seqs: set[int], cluster_group: set[Cluster],
        seq_assignments: list[list[Cluster]], l: int,
        nkl: np.ndarray = None, emission: int = None,
    ):
        self.uuid = uuid.uuid4().int

        self.seqs = seqs
        self.cluster_group = cluster_group
        self.cluster_group.add(self)
        self.seq_assignments = seq_assignments
        for i in seqs:
            seq_assignments[i][l] = self
        self.l = l

        self.emission = emission
        self.nkl = nkl
        if not self.emission is None:
            self.nkl[emission] += 1

        self.children = set()
        self.parents = set()

    def __hash__(self) -> int:
        return self.uuid

    def __eq__(self, other) -> bool:
        return isinstance(other, Cluster) and (self.uuid == other.uuid)

    def __len__(self) -> int:
        return len(self.seqs)

    def __repr__(self) -> str:
        return str(self.seqs)

    def add_child(self, child: Cluster):
        self.children.add(child)
        child.parents.add(self)

    def remove(self, idx: int):
        self.seqs.remove(idx)
        self.seq_assignments[idx][self.l] = None
        if len(self.seqs) == 0:
            for parent in self.parents:
                parent.children.remove(self)
            for child in self.children:
                child.parents.remove(self)
            if not self.emission is None:
                self.nkl[self.emission] -= 1
            self.cluster_group.remove(self)

    def add(self, idx: int):
        self.seqs.add(idx)
        assert self.seq_assignments[idx][self.l] is None
        self.seq_assignments[idx][self.l] = self


@profile
def me(x: np.ndarray, HP: HyperParams, viz: bool = False):
    # Init parameters.
    mu_alpha = HP.tau_1 / HP.tau_2
    sigma2_alpha = mu_alpha / HP.tau_2
    mu_log_alpha = special.digamma(HP.tau_1) - np.log(HP.tau_2)
    sigma2_log_alpha = special.polygamma(1, HP.tau_1)

    mu_d = np.full(HP.L-1, HP.v_1 / (HP.v_1 + HP.v_2))
    sigma2_d = np.full(HP.L-1, (HP.v_1 * HP.v_2) / ((HP.v_1 + HP.v_2)**2 * (HP.v_1 + HP.v_2 + 1)))
    mu_log_d = np.full(HP.L-1, special.digamma(HP.v_1) - special.digamma(HP.v_1 + HP.v_2))
    sigma2_logit_d = np.full(HP.L-1, float("-inf"))

    mu_gamma = np.full(HP.L, HP.phi_1 / HP.phi_2)
    sigma2_gamma = mu_gamma / HP.phi_2
    mu_log_gamma = np.full(HP.L, special.digamma(HP.phi_1) - np.log(HP.phi_2))
    sigma2_log_gamma = np.full(HP.L, special.polygamma(1, HP.phi_1))

    # Init clusters.
    nk = np.zeros((HP.L, HP.K), dtype=np.int64)
    r_assignments = [[None for _ in range(HP.L)] for _ in range(HP.N)]
    q_assignments = [[None for _ in range(HP.L-1)] for _ in range(HP.N)]
    rs = [set() for _ in range(HP.L)]
    qs = [set() for _ in range(HP.L-1)]
    x_modes = stats.mode(x, axis=0).mode
    r = Cluster(seqs=set(range(HP.N)), cluster_group=rs[0],
                seq_assignments=r_assignments, l=0,
                nkl=nk[0], emission=x_modes[0])
    for l in range(HP.L-1):
        q = Cluster(seqs=set(range(HP.N)), cluster_group=qs[l],
                    seq_assignments=q_assignments, l=l)
        r.add_child(q)
        r = Cluster(seqs=set(range(HP.N)), cluster_group=rs[l+1],
                    seq_assignments=r_assignments, l=l+1,
                    nkl=nk[l+1], emission=x_modes[l+1])
        q.add_child(r)

    # ME.
    if viz: log = Logger()
    early_stop = EarlyStopping(patience=5, minimize=False)
    while not early_stop.converged():
        max_step(
            x, HP,
            mu_alpha, sigma2_alpha, mu_log_alpha,
            mu_d, sigma2_d, mu_log_d,
            mu_gamma, sigma2_gamma,
            rs, qs, nk, r_assignments, q_assignments
        )
        mu_alpha, sigma2_alpha, mu_log_alpha, sigma2_log_alpha = expect_step(
            HP,
            mu_d, sigma2_d, mu_log_d, sigma2_logit_d,
            mu_gamma, sigma2_gamma, mu_log_gamma, sigma2_log_gamma,
            rs, qs, nk
        )
        elbo = calc_elbo(
            HP,
            mu_alpha, sigma2_alpha, mu_log_alpha, sigma2_log_alpha,
            mu_d, sigma2_d, mu_log_d, sigma2_logit_d,
            mu_gamma, sigma2_gamma, mu_log_gamma, sigma2_log_gamma,
            rs, qs, nk
        )
        print(f"{early_stop.step}.  elbo={elbo:.4f}")
        early_stop.update(elbo)
        if viz:
            log.log({
                "elbo": elbo,
                "mu_alpha": mu_alpha, "sigma2_alpha": sigma2_alpha,
                "mu_d": mu_d, "sigma2_d": sigma2_d,
                "mu_gamma": mu_gamma, "sigma2_gamma": sigma2_gamma,
            })

    if viz:
        draw_viz(HP, rs, qs, "clusters")
        log.plot("log.png")


@profile
def calc_elbo(
    HP: HyperParams,
    mu_alpha: float, sigma2_alpha: float, mu_log_alpha: float, sigma2_log_alpha: float,
    mu_d: np.ndarray, sigma2_d: np.ndarray, mu_log_d: np.ndarray, sigma2_logit_d: np.ndarray,
    mu_gamma: np.ndarray, sigma2_gamma: np.ndarray, mu_log_gamma: np.ndarray, sigma2_log_gamma: np.ndarray,
    rs: list[set[Cluster]], qs: list[set[Cluster]],
    nk: np.ndarray,
) -> float:
    elbo = delta_ElogGamma(mu_alpha, sigma2_alpha) - delta_ElogGamma(mu_alpha, sigma2_alpha, b=HP.N)
    elbo += (sum(len(r) for r in rs)+HP.tau_1)*mu_log_alpha - HP.tau_2*mu_alpha  # -1 cancels out w/ alpha entropy.
    elbo += HP.tau_1*np.log(HP.tau_2) - special.gammaln(HP.tau_1)

    for l in range(HP.L-1):
        # -1s cancel out in dl entropy.
        elbo += (len(qs[l])-len(rs[l])-len(rs[l+1])+HP.v_1) * mu_log_d[l]
        elbo += HP.v_2 * delta_Elogx(mu_d[l], sigma2_d[l], a=-1, b=1)
        elbo -= special.betaln(HP.v_1, HP.v_2)

        elbo -= len(qs[l]) * delta_ElogGamma(mu_d[l], sigma2_d[l], a=-1, b=1)
        elbo += sum(delta_ElogGamma(mu_d[l], sigma2_d[l], a=-1, b=len(b)) for b in qs[l])

        # alpha and d term.
        z = mu_alpha / mu_d[l]
        dd2 = special.digamma(z)*(2*z/mu_d[l]**2) + special.polygamma(1, z)*z**2/mu_d[l]**2
        elbo += delta_ElogGamma(mu_alpha, sigma2_alpha, a=1/mu_d[l]) + 0.5*sigma2_d[l]*dd2

        z += len(qs[l])
        dd2 = special.digamma(z)*(2*mu_alpha/mu_d[l]**3) + special.polygamma(1, z)*mu_alpha**2/mu_d[l]**4
        elbo -= delta_ElogGamma(mu_alpha, sigma2_alpha, a=1/mu_d[l], b=len(qs[l])) + 0.5*sigma2_d[l]*dd2

    elbo += (HP.phi_1*mu_log_gamma - HP.phi_2*mu_gamma).sum()  # -1 cancels out w/ gammal entropy.
    for l in range(HP.L):
        elbo += HP.phi_1*np.log(HP.phi_2) - special.gammaln(HP.phi_1)
        elbo += delta_ElogGamma(mu_gamma[l], sigma2_gamma[l], a=HP.K)
        elbo -= delta_ElogGamma(mu_gamma[l], sigma2_gamma[l], a=HP.K, b=len(rs[l]))
        elbo += sum(delta_ElogGamma(mu_gamma[l], sigma2_gamma[l], b=n) for n in nk[l])
        elbo -= HP.K*delta_ElogGamma(mu_gamma[l], sigma2_gamma[l])

    elbo += sum(special.gammaln(len(a)) for a in rs[0])
    for l in range(HP.L-1):
        elbo += sum(special.gammaln(len(a.children)) - special.gammaln(len(a)) for a in rs[l])
        elbo += sum(special.gammaln(len(a.parents)) for a in rs[l+1])

    # Variational entropy term.
    elbo += normal_entropy(sigma2_log_alpha)
    elbo += normal_entropy(sigma2_log_gamma).sum() + normal_entropy(sigma2_logit_d).sum()
    return elbo


def delta_ElogGamma(mu: float, sigma2: float, a: float = 1.0, b: float = 0.0) -> float:
    x = a*mu + b
    return special.gammaln(x) + 0.5*sigma2*a**2 * special.polygamma(1, x)


def normal_entropy(sigma2: float) -> float:
    return 0.5 * np.log(2*np.pi*np.e * sigma2)


class Logger:
    def __init__(self):
        self.vals = {}

    def log(self, log_dict: dict) -> None:
        for (k, v) in log_dict.items():
            if k not in self.vals:
                self.vals[k] = []
            self.vals[k].append(v.copy())

    def plot(self, save_path: str) -> None:
        fig, axs = plt.subplots(ncols=len(self.vals), figsize=(len(self.vals)*8, 6))
        fig.tight_layout(pad=2)
        for (k, v), ax in zip(self.vals.items(), axs):
            if isinstance(v[0], float):
                ax.plot(v, color="blue")
            elif isinstance(v[0], np.ndarray):
                for i, vi in enumerate(v):
                    ax.plot(vi, alpha=0.4 + 0.6*(i/(len(v)-1)), color="blue")
            else:
                raise ValueError(f"Logger does not support {k} of type {type(v)}")
            ax.set_title(k)
        fig.savefig(save_path)


class EarlyStopping:
    def __init__(self, patience: int, minimize: bool = True, tol: float = 1e-5):
        self.patience = patience
        self.minimize = minimize
        self.tol = tol

        self.step = 0

        self.min_val = float("inf")
        self.steps_since_min = 0

    def update(self, x: float):
        self.step += 1
        if not self.minimize:
            x = -x
        if self.min_val - x > self.tol:
            self.min_val = x
            self.steps_since_min = 0
        else:
            self.steps_since_min += 1

    def converged(self) -> bool:
        return self.steps_since_min > self.patience


def draw_viz(HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]], save_name: str) -> None:
    d = graphviz.Digraph("dfcp", graph_attr={"rankdir": "LR"})
    cluster_names = {}
    for l in reversed(range(HP.L)):
        for i, r in enumerate(rs[l]):
            name = f"R,{l},{i}"
            text = f"{name}\nn={len(r.seqs)}\nx={r.emission}"
            d.node(name, text, style="filled", fillcolor="lightblue")
            cluster_names[r] = name
        if l == HP.L-1:
            continue

        for i, q in enumerate(qs[l]):
            name = f"Q,{l},{i}"
            text = f"{name}\nn={len(q.seqs)}"
            d.node(name, text, style="filled", fillcolor="lightgreen", shape="diamond")
            cluster_names[q] = name
            d.edges([(name, cluster_names[child]) for child in q.children])

        for r in rs[l]:
            d.edges([(cluster_names[r], cluster_names[child]) for child in r.children])
    d.render(save_name, format="png", cleanup=True)


@profile
def max_step(
    x: np.ndarray, HP: HyperParams,
    mu_alpha: float, sigma2_alpha: float, mu_log_alpha: float,
    mu_d: np.ndarray, sigma2_d: np.ndarray, mu_log_d: np.ndarray,
    mu_gamma: np.ndarray, sigma2_gamma: np.ndarray,
    rs: list[set[Cluster]], qs: list[set[Cluster]],
    nk: np.ndarray,
    r_assignments: list[list[Cluster]], q_assignments: list[list[Cluster]],
) -> None:
    for i in range(HP.N):
        for a in r_assignments[i]:
            a.remove(i)
        for b in q_assignments[i]:
            b.remove(i)

        messages = []
        for l in reversed(range(HP.L)):
            # Compute likelihood term for a-message.
            ma = {}
            log_likelihood = (
                delta_Elogx(mu_gamma[l], sigma2_gamma[l], b=nk[l, x[i, l]])
                - delta_Elogx(mu_gamma[l], sigma2_gamma[l], a=HP.K, b=len(rs[l]))
            )
            ma["new"] = (log_likelihood, None)
            for a in rs[l]:
                ma[a] = (0 if x[i, l] == a.emission else float("-infinity"), None)
            if l == HP.L - 1:
                messages.append((ma,))
                continue

            # b-messages.
            next_ma = messages[-1][0]
            mb = {}
            for b in qs[l]:
                assert len(b.children) == 1
                next_a = next(iter(b.children))
                mb[b] = (next_ma[next_a][0], next_a)

            nQl = len(qs[l])
            mu_y = mu_alpha + nQl*mu_d[l]
            sigma2_y = sigma2_alpha + nQl*nQl * sigma2_d[l]
            elogy = delta_Elogx(mu_y, sigma2_y)

            best_a = "new"
            best_a_ll = mu_log_alpha + next_ma["new"][0]
            for a in rs[l+1]:
                ll = mu_log_d[l] + np.log(len(a.parents)) + next_ma[a][0]
                if ll > best_a_ll:
                    best_a = a
                    best_a_ll = ll
            mb["new"] = (-elogy + best_a_ll, best_a)

            # a-messages.
            ma["new"] = (ma["new"][0] + mb["new"][0], "new")
            for a in rs[l]:
                best_b = "new"
                best_b_ll = np.log(len(a.children)) + mu_log_d[l] + mb["new"][0]
                for b in a.children:
                    ll = delta_Elogx(mu_d[l], sigma2_d[l], a=-1, b=len(b)) + mb[b][0]
                    if ll > best_b_ll:
                        best_b = b
                        best_b_ll = ll
                ma[a] = (ma[a][0] - np.log(len(a)) + best_b_ll, best_b)

            messages.append((ma, mb))

        # Backtrack for viterbi path.
        path = []
        for l, (ma, mb) in enumerate(reversed(messages[1:])):
            if l == 0:
                a = max(ma, key=lambda x: ma.get(x)[0])
                path.append(a)

            b = ma[a][1]
            path.append(b)

            a = mb[b][1]
            path.append(a)

        for j in range(len(path)):
            l = j // 2
            if path[j] != "new":
                path[j].add(i)
                continue

            if j % 2 == 0:
                new_cluster = Cluster(seqs=set([i]), cluster_group=rs[l],
                                      seq_assignments=r_assignments, l=l,
                                      nkl=nk[l], emission=x[i, l])
            else:
                new_cluster = Cluster(seqs=set([i]), cluster_group=qs[l],
                                      seq_assignments=q_assignments, l=l)
            if j - 1 >= 0:
                path[j-1].add_child(new_cluster)
            if j + 1 < len(path) and path[j+1] != "new":
                new_cluster.add_child(path[j+1])
            path[j] = new_cluster


@profile
def expect_step(
    HP: HyperParams,
    mu_d: np.ndarray, sigma2_d: np.ndarray, mu_log_d: np.ndarray, sigma2_logit_d: np.ndarray,
    mu_gamma: np.ndarray, sigma2_gamma: np.ndarray, mu_log_gamma: np.ndarray, sigma2_log_gamma: np.ndarray,
    rs: list[set[Cluster]], qs: list[set[Cluster]],
    nk: np.ndarray
) -> tuple[float, float, float, float]:
    # Expectation.
    # alpha update.
    mu_alpha, sigma2_alpha, mu_log_alpha, sigma2_log_alpha = laplace_log_approx(
        ll_alpha, ll_log_alpha_d2, args=(mu_d, sigma2_d, HP, rs, qs)
    )

    for l in range(HP.L):
        # gamma update.
        mu_gamma[l], sigma2_gamma[l], mu_log_gamma[l], sigma2_log_gamma[l] = laplace_log_approx(
            ll_gammal, ll_log_gammal_d2, args=(HP, len(rs[l]), nk[l])
        )

        if l >= HP.L-1:
            break
        # d update.
        res = optimize.minimize_scalar(
            lambda eta: -ll_logit_dl(eta, mu_alpha, sigma2_alpha, HP, rs, qs, l),
            method="bounded", bounds=(-10, 10)
        )
        assert res.success, res.message
        eta_mode = res.x
        d_mode = special.expit(eta_mode)
        eta_var = -1/ll_logit_dl_d2(d_mode, mu_alpha, sigma2_alpha, HP, rs, qs, l)
        assert eta_var > 0, eta_var

        mu_d[l] = d_mode + 0.5 * eta_var * (1-2*d_mode) * d_mode * (1-d_mode)
        sigma2_d[l] = (
            d_mode**2
            + 0.5 * eta_var * (4*d_mode - 6*d_mode**2) * d_mode * (1-d_mode)
            - mu_d[l]**2
        )
        assert sigma2_d[l] > 0, sigma2_d[l]
        mu_log_d[l] = delta_Elogx(mu_d[l], sigma2_d[l])
        sigma2_logit_d[l] = eta_var
    return mu_alpha, sigma2_alpha, mu_log_alpha, sigma2_log_alpha


def delta_Elogx(mu: float, sigma2: float, a: float = 1.0, b: float = 0.0) -> float:
    x = a*mu + b
    return np.log(x) - 0.5 * sigma2 * a**2 / x**2


def laplace_log_approx(ll_func, ll_log_d2_func, args: tuple = ()) -> tuple[float, float, float, float]:
    res = optimize.minimize_scalar(
        lambda eta: -(ll_func(np.exp(eta), *args) + eta),
        method="bounded", bounds=(-10, 10)
    )
    assert res.success, res.message
    eta_mode = res.x

    eta_var = -1/ll_log_d2_func(np.exp(eta_mode), *args)
    assert eta_var > 0, eta_var

    mu = np.exp(eta_mode + eta_var/2)
    sigma2 = (np.exp(eta_var) - 1) * np.exp(2*eta_mode + eta_var)
    assert sigma2 > 0, sigma2
    return mu, sigma2, eta_mode, eta_var


def ll_alpha(
    alpha: float, mu_d: np.ndarray, sigma2_d: np.ndarray,
    HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]]
) -> float:
    ll = special.gammaln(alpha) - special.gammaln(alpha + HP.N)
    ll += (sum([len(r) for r in rs]) + HP.tau_1 - 1) * np.log(alpha) 
    ll -= HP.tau_2*alpha
    for l in range(HP.L-1):
        idxs = np.arange(len(qs[l]))
        ll -= (
            np.log(alpha/mu_d[l] + idxs)
            + 0.5*sigma2_d[l] * (alpha**2 + 2*alpha*idxs*mu_d[l]) / (mu_d[l]**2 * (alpha + idxs*mu_d[l])**2)
        ).sum()
    return ll


def ll_log_alpha_d2(
    alpha: float, mu_d: np.ndarray, sigma2_d: np.ndarray,
    HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]]
) -> float:
    d2 = (1 / (alpha + np.arange(HP.N))**2).sum()
    for l in range(HP.L-1):
        idxs = np.arange(len(qs[l]))
        z = alpha + idxs*mu_d[l]
        d2 += (1/z**2 + 3*sigma2_d[l]*idxs**2 / z**4).sum()
    d2 *= alpha**2
    d2 -= HP.tau_1 + sum([len(r) for r in rs])
    return d2


def ll_gammal(gamma: float, HP: HyperParams, nRl: int, nkl: np.ndarray) -> float:
    ll = (HP.phi_1-1)*np.log(gamma) - HP.phi_2*gamma
    ll += special.gammaln(HP.K*gamma) - special.gammaln(HP.K*gamma + nRl)
    ll += special.gammaln(gamma+nkl).sum() - HP.K*special.gammaln(gamma)
    return ll


def ll_log_gammal_d2(gamma: float, HP: HyperParams, nRl: int, nkl: np.ndarray) -> float:
    d2 = -HP.phi_1 + gamma*gamma * (
        (HP.K**2 / (HP.K*gamma + np.arange(nRl))**2).sum()
        - sum([(1/(gamma + np.arange(n))**2).sum() for n in nkl])
    )
    return d2


def ll_logit_dl(
    eta: float, mu_alpha: float, sigma2_alpha: float,
    HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]], l: int
) -> float:
    d = special.expit(eta)
    return ll_dl(d, mu_alpha, sigma2_alpha, HP, rs, qs, l) + np.log(d) + np.log(1-d)


def ll_dl(
    d: float, mu_alpha: float, sigma2_alpha: float,
    HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]], l: int
) -> float:
    nQl = len(qs[l])
    ll = (nQl - len(rs[l]) - len(rs[l+1]) + HP.v_1 - 1) * np.log(d)
    ll += (HP.v_2 - 1) * np.log(1 - d)
    ll += -nQl * special.loggamma(1-d) + sum([special.loggamma(len(b) - d) for b in qs[l]])
    idxs = np.arange(nQl)
    ll -= (np.log(mu_alpha/d + idxs) - sigma2_alpha/(2* (mu_alpha+idxs*d)**2)).sum()
    return ll


def ll_logit_dl_d2(
    d: float, mu_alpha: float, sigma2_alpha: float,
    HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]], l: int
) -> float:
    nQl = len(qs[l])
    d2 = (len(rs[l]) - nQl + len(rs[l+1]) + 1 - HP.v_1) / (d*d)
    d2 -= (HP.v_2 - 1) / (1 - d)**2
    d2 += -nQl * special.polygamma(1, 1-d) + sum([special.polygamma(1, len(b) - d) for b in qs[l]])
    idxs = np.arange(nQl)
    z = mu_alpha + idxs*d
    d2 -= (
        (mu_alpha*(mu_alpha+2*idxs*d)) / (d*d* z**2)
        - (3*idxs*idxs*sigma2_alpha / z**4)
    ).sum()

    d2 *= (d*(1-d)) ** 2
    d2 -= (1-2*d)**2 + 2*d*(1-d)
    return d2


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--seed", default=0, type=int)
    parser.add_argument("--viz", action="store_true")
    HP = HyperParams().add_params(parser)
    args = parser.parse_args()
    HP.override_args(args)

    x, _, _ = generate(HP, seed=args.seed)
    print(f"x:\n{x}")

    me(x, HP, viz=args.viz)

