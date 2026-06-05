from __future__ import annotations
import uuid
from pprint import pprint

import numpy as np
from scipy import stats, special, optimize
from line_profiler import profile

from generate import generate, HyperParams


class Cluster:
    def __init__(
        self, seqs: set[int], cluster_group: set[Cluster],
        seq_assignments: list[list[Cluster]], l: int,
        seg_group: set[Cluster] = None, emission: int = None,
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
        self.seg_group = seg_group
        if not self.emission is None:
            self.seg_group[emission].add(self)

        self.children = set()
        self.parents = set()

    def __hash__(self) -> int:
        return self.uuid

    def __eq__(self, other) -> bool:
        return isinstance(other, Cluster) and (self.uuid == other.uuid)

    def __len__(self) -> int:
        return len(self.seqs)

    def __repr__(self) -> str:
        return f"E {self.emission}   {self.seqs}"

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
                self.seg_group[self.emission].remove(self)
            self.cluster_group.remove(self)

    def add(self, idx: int):
        self.seqs.add(idx)
        assert self.seq_assignments[idx][self.l] is None
        self.seq_assignments[idx][self.l] = self


@profile
def me(x: np.ndarray, HP: HyperParams):
    # Init parameters.
    mu_alpha = HP.tau_1 / HP.tau_2
    sigma2_alpha = mu_alpha / HP.tau_2

    mu_d = np.full(HP.L-1, HP.v_1 / (HP.v_1 + HP.v_2))
    sigma2_d = np.full(HP.L-1, (HP.v_1 * HP.v_2) / ((HP.v_1 + HP.v_2)**2 * (HP.v_1 + HP.v_2 + 1)))

    mu_gamma = np.full(HP.L, HP.phi_1 / HP.phi_2)
    sigma2_gamma = mu_gamma / HP.phi_2

    # Init clusters.
    segregated_rs = [[set() for _ in range(HP.K)] for _ in range(HP.L)]
    r_assignments = [[None for _ in range(HP.L)] for _ in range(HP.N)]
    q_assignments = [[None for _ in range(HP.L-1)] for _ in range(HP.N)]
    rs = [set() for _ in range(HP.L)]
    qs = [set() for _ in range(HP.L-1)]
    x_modes = stats.mode(x, axis=0).mode
    r = Cluster(seqs=set(range(HP.N)), cluster_group=rs[0],
                seq_assignments=r_assignments, l=0,
                seg_group=segregated_rs[0], emission=x_modes[0])
    for l in range(HP.L-1):
        q = Cluster(seqs=set(range(HP.N)), cluster_group=qs[l],
                    seq_assignments=q_assignments, l=l)
        r.add_child(q)
        r = Cluster(seqs=set(range(HP.N)), cluster_group=rs[l+1],
                    seq_assignments=r_assignments, l=l+1,
                    seg_group=segregated_rs[l+1], emission=x_modes[l+1])
        q.add_child(r)
    print("rs:")
    pprint(rs)
    print("qs:")
    pprint(qs)

    # Maximization.
    for i in range(HP.N):
        for a in r_assignments[i]:
            a.remove(i)
        for b in q_assignments[i]:
            b.remove(i)

        messages = []
        for l in reversed(range(HP.L)):
            # Compute likelihood term for a-message.
            ma = {}
            nkl = len(segregated_rs[l][x[i, l]])
            nRl = len(rs[l])
            log_likelihood = (
                delta_Elogx(mu_gamma[l], sigma2_gamma[l], b=nkl)
                - delta_Elogx(mu_gamma[l], sigma2_gamma[l], a=HP.K, b=nRl)
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

            next_a_ll = {}
            # TODO: EE[log a] is closed form for log normal.
            next_a_ll["new"] = delta_Elogx(mu_alpha, sigma2_alpha) + next_ma["new"][0]
            for a in rs[l+1]:
                next_a_ll[a] = delta_Elogx(mu_d[l], sigma2_d[l]) + np.log(len(a.parents)) + next_ma[a][0]
            next_a = max(next_a_ll, key=next_a_ll.get)
            mb["new"] = (-elogy + next_a_ll[next_a], next_a)

            # a-messages.
            ma["new"] = (ma["new"][0] + mb["new"][0], "new")
            for a in rs[l]:
                next_b_ll = {}
                next_b_ll["new"] = np.log(len(a.children)) + delta_Elogx(mu_d[l], sigma2_d[l]) + mb["new"][0]
                for b in a.children:
                    next_b_ll[b] = delta_Elogx(mu_d[l], sigma2_d[l], a=-1, b=len(b)) + mb[b][0]
                next_b = max(next_b_ll, key=next_b_ll.get)
                ma[a] = (ma[a][0] - np.log(len(a)) + next_b_ll[next_b], next_b)

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
                                      seg_group=segregated_rs[l], emission=x[i, l])
            else:
                new_cluster = Cluster(seqs=set([i]), cluster_group=qs[l],
                                      seq_assignments=q_assignments, l=l)
            if j - 1 >= 0:
                path[j-1].add_child(new_cluster)
            if j + 1 < len(path) and path[j+1] != "new":
                new_cluster.add_child(path[j+1])
            path[j] = new_cluster
    print("rs:")
    pprint(rs)
    print("qs:")
    pprint(qs)

    # Expectation.
    # alpha update.
    mu_alpha, sigma2_alpha = laplace_log_approx(ll_alpha, ll_log_alpha_d2,
                                                args=(mu_d, sigma2_d, HP, rs, qs))
    print(f"mu_alpha: {mu_alpha}")
    print(f"sigma2_alpha: {sigma2_alpha}")

    for l in range(HP.L):
        # gamma update.
        mu_gamma[l], sigma2_gamma[l] = laplace_log_approx(ll_gammal, ll_log_gammal_d2,
                                                          args=(HP, len(rs[l]), segregated_rs[l]))

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

    print(f"mu_gamma: {mu_gamma}")
    print(f"sigma2_gamma: {sigma2_gamma}")
    print(f"mu_d: {mu_d}")
    print(f"sigma2_d: {sigma2_d}")


def delta_Elogx(mu: float, sigma2: float, a: float = 1.0, b: float = 0.0) -> float:
    x = a*mu + b
    return np.log(x) - 0.5 * sigma2 * a**2 / x**2


def laplace_log_approx(ll_func, ll_log_d2_func, args: tuple = ()) -> tuple[float, float]:
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
    return mu, sigma2


def ll_alpha(
    alpha: float, mu_d: np.ndarray, sigma2_d: np.ndarray,
    HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]]
) -> float:
    ll = -np.log(alpha + np.arange(HP.N)).sum()
    ll += sum([len(r) * np.log(alpha) for r in rs])
    for l in range(HP.L-1):
        idxs = np.arange(len(qs[l]))
        ll -= (
            np.log(alpha/mu_d[l] + idxs)
            + 0.5*sigma2_d[l] * (alpha**2 + 2*alpha*idxs*mu_d[l]) / (mu_d[l]**2 * (alpha + idxs*mu_d[l])**2)
        ).sum()
    ll += (HP.tau_1-1)*np.log(alpha) - HP.tau_2*alpha
    return ll


def ll_log_alpha_d2(
    alpha: float, mu_d: np.ndarray, sigma2_d: np.ndarray,
    HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]]
) -> float:
    d2 = (alpha**2 / (alpha + np.arange(HP.N))**2).sum()
    d2 -= sum([len(r) for r in rs])
    for l in range(HP.L-1):
        idxs = np.arange(len(qs[l]))
        d2 += alpha**2 * (
            1 / (alpha + idxs*mu_d[l])**2
            + 3*sigma2_d[l]*idxs**2 / (alpha+idxs*mu_d[l])**4
        ).sum()
    d2 -= HP.tau_1
    return d2


# def d2_finite_diff(f, x: float, eps: float = 1e-9) -> float:
#     return (f(x+2*eps) - 2*f(x+eps) + f(x)) / (eps*eps)


def ll_gammal(gamma: float, HP: HyperParams, nRl: int, segregated_l: list[set[Cluster]]) -> float:
    ll = (HP.phi_1-1)*np.log(gamma) - HP.phi_2*gamma
    ll -= np.log(HP.K*gamma + np.arange(nRl)).sum()
    ll += sum([np.log(gamma + np.arange(len(seg))).sum() for seg in segregated_l])
    return ll


def ll_log_gammal_d2(gamma: float, HP: HyperParams, nRl: int, segregated_l: list[set[Cluster]]) -> float:
    d2 = -HP.phi_1 + gamma*gamma * (
        (HP.K**2 / (HP.K*gamma + np.arange(nRl))**2).sum()
        - sum([(1/(gamma + np.arange(len(seg)))**2).sum() for seg in segregated_l])
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
    d2 -= (
        (mu_alpha*(mu_alpha+2*idxs*d))/(d*d* (mu_alpha + idxs*d)**2)
        - (3*idxs*idxs*sigma2_alpha/(mu_alpha+idxs*d)**4)
    ).sum()

    d2 *= (d*(1-d)) ** 2
    d2 -= (1-2*d)**2 + 2*d*(1-d)
    return d2


if __name__ == "__main__":
    # TODO: different HP for generate and me.
    HP = HyperParams().cli()
    x, _, _ = generate(HP)
    print(f"x:\n{x}")

    me(x, HP)

