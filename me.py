from __future__ import annotations
import uuid
from pprint import pprint

import numpy as np
from scipy.stats import gamma, beta, mode, lognorm
from scipy.optimize import minimize_scalar

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

    def __hash__(self):
        return self.uuid

    def __eq__(self, other):
        return isinstance(other, Cluster) and (self.uuid == other.uuid)

    def __repr__(self):
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
            if self.emission:
                self.seg_group[self.emission].remove(self)
            self.cluster_group.remove(self)

    def add(self, idx: int):
        self.seqs.add(idx)
        assert self.seq_assignments[idx][self.l] is None
        self.seq_assignments[idx][self.l] = self


def me(x: np.ndarray, HP: HyperParams):
    # Init parameters.
    alpha = gamma(a=HP.tau_1, scale=1/HP.tau_2)
    d = [beta(a=HP.v_1, b=HP.v_2) for _ in range(HP.L-1)]
    gammal = [gamma(a=HP.phi_1, scale=1/HP.phi_2) for _ in range(HP.L)]

    # Init clusters.
    segregated_rs = [[set() for _ in range(HP.K)] for _ in range(HP.L)]
    r_assignments = [[None for _ in range(HP.L)] for _ in range(HP.N)]
    q_assignments = [[None for _ in range(HP.L-1)] for _ in range(HP.N)]
    rs = [set() for _ in range(HP.L)]
    qs = [set() for _ in range(HP.L-1)]
    x_modes = mode(x, axis=0)[0]
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

            ma["new"] = (gammal[l].expect(lambda g: np.log(g + nkl) - np.log(HP.K * g + nRl)), None)
            for a in rs[l]:
                ma[a] = (0 if x[i, l] == a.emission else float("-infinity"), None)
            if l == HP.L - 1:
                messages.append((ma,))
                continue

            mb = {}
            for b in qs[l]:
                assert len(b.children) == 1
                next_a = next(iter(b.children))
                mb[b] = (messages[-1][0][next_a][0], next_a)

            nQl = len(qs[l])
            mu_y = alpha.mean() + nQl * d[l].mean()
            sigma2_y = alpha.var() + nQl*nQl * d[l].var()
            elogy = np.log(mu_y) - 0.5 * sigma2_y / (mu_y*mu_y)
            next_ma = {}
            next_ma["new"] = alpha.expect(np.log) + messages[-1][0]["new"][0]
            for a in rs[l+1]:
                next_ma[a] = d[l].expect(np.log) + np.log(len(a.parents)) + messages[-1][0][a][0]
            next_a = max(next_ma, key=next_ma.get)
            mb["new"] = (-elogy + next_ma[next_a], next_a)

            ma["new"] = (ma["new"][0] + mb["new"][0], "new")
            for a in rs[l]:
                next_mb = {}
                next_mb["new"] = np.log(len(a.children)) + d[l].expect(np.log) + mb["new"][0]
                for b in a.children:
                    next_mb[b] = d[l].expect(lambda x: np.log(len(b.seqs) - x)) + mb[b][0]
                next_b = max(next_mb, key=next_mb.get)
                ma[a] = (ma[a][0] - np.log(len(a.seqs)) + next_mb[next_b], next_b)

            messages.append((ma, mb))

        path = []
        for l, (ma, mb) in enumerate(reversed(messages[1:])):
            if l == 0:
                a = max(ma, key=lambda x: ma.get(x)[0])
                path.append(a)

            if l == HP.L-1:
                break
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
    res = minimize_scalar(lambda eta: -(ll_alpha(np.exp(eta), d, HP, rs, qs) + eta),
                          method="bounded", bounds=(-10, 10))
    assert res.success, res.message
    eta_mode = res.x
    print(eta_mode)
    eta_var = -1/ll_log_alpha_d2(np.exp(eta_mode), d, HP, rs, qs)
    assert eta_var > 0, eta_var
    print(eta_var)
    alpha = lognorm(s=np.sqrt(eta_var), scale=np.exp(eta_mode))

    # TODO: gamma updates are wrong. eta_var < 0.
    for l in range(HP.L):
        res = minimize_scalar(lambda eta: -(ll_gammal(np.exp(eta), HP, segregated_rs[l]) + eta),
                              method="bounded", bounds=(-10, 10))
        assert res.success, res.message
        eta_mode = res.x
        print(eta_mode)
        eta_var = -1/ll_log_gammal_d2(np.exp(eta_mode), HP, segregated_rs[l])
        assert eta_var > 0, eta_var
        print(eta_var)
        gammal[l] = lognorm(s=np.sqrt(eta_var), scale=np.exp(eta_mode))


def ll_alpha(
    alpha: float, d: np.ndarray, HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]]
) -> float:
    ll = -np.log(alpha + np.arange(HP.N)).sum()
    ll += sum([len(r) * np.log(alpha) for r in rs])
    for l in range(1, HP.L-1):
        idxs = np.arange(len(qs[l]))
        ll -= (
            np.log(alpha/d[l].mean()+idxs)
            +0.5*d[l].var()*alpha*alpha+2*alpha*idxs*d[l].mean()
            /(d[l].mean()**2 * (alpha+idxs*d[l].mean())**2)
        ).sum()
    ll += (HP.tau_1-1)*np.log(alpha) - HP.tau_2*alpha
    return ll


def ll_log_alpha_d2(
    alpha: float, d: np.ndarray, HP: HyperParams, rs: list[set[Cluster]], qs: list[set[Cluster]]
) -> float:
    d2 = (alpha*alpha / (alpha + np.arange(HP.N))**2).sum()
    d2 -= sum([len(r) for r in rs])
    for l in range(1, HP.L-1):
        idxs = np.arange(len(qs[l]))
        d2 += alpha*alpha * (
            1/(alpha+idxs*d[l].mean())**2
            + 2*d[l].var()*idxs
            / (d[l].mean()*(alpha+idxs*d[l].mean())**3)
        ).sum()
    d2 -= HP.tau_1
    return d2


# def d2_finite_diff(f, x: float, eps: float = 1e-9) -> float:
#     return (f(x+2*eps) - 2*f(x+eps) + f(x)) / (eps*eps)


def ll_gammal(gamma: float, HP: HyperParams, segregated_l: list[set[Cluster]]) -> float:
    ll = (HP.phi_1-1)*np.log(gamma) - HP.phi_2*gamma
    ll -= np.log(4*gamma + np.arange(HP.N)).sum()
    ll += sum([np.log(gamma + np.arange(len(seg))).sum() for seg in segregated_l])
    return ll


def ll_log_gammal_d2(gamma: float, HP: HyperParams, segregated_l: list[set[Cluster]]) -> float:
    d2 = HP.tau_1 + gamma*gamma * (
        (16 / (4*gamma + np.arange(HP.N))).sum()
        - sum([((gamma + np.arange(len(seg)))**-2).sum() for seg in segregated_l])
    )
    return d2


if __name__ == "__main__":
    # TODO: different HP for generate and me.
    HP = HyperParams().cli()
    x, _, _ = generate(HP)
    print(f"x:\n{x}")

    me(x, HP)
