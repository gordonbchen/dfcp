from __future__ import annotations
import random
from argparse import ArgumentParser, Namespace
from dataclasses import dataclass, asdict

import numpy as np
from scipy.stats import gamma, beta, dirichlet


class CLIParams:
    def add_params(self, parser: ArgumentParser) -> CLIParams:
        for k, v in asdict(self).items():
            parser.add_argument(f"--{k}", type=type(v), default=v)
        return self

    def override_args(self, args: Namespace) -> CLIParams:
        for k in asdict(self):
            setattr(self, k, vars(args)[k])
        return self


@dataclass
class HyperParams(CLIParams):
    N: int = 8   # Num sequences.
    L: int = 16  # Seq length.
    K: int = 4   # Num alleles.

    # alpha Gamma hyperparams (alpha, rate).
    tau_1: float = 1.0
    tau_2: float = 1.0

    # d Beta hyperparams.
    v_1: float = 1.0
    v_2: float = 1.0

    # gamma Gamma hyperparams.
    phi_1: float = 1.0
    phi_2: float = 1.0


def crp(nums: list[int], alpha: float, d: float) -> list[list[int]]:
    clusters = [[nums[0]]]
    discounted_sizes = [1-d]
    for i, num in enumerate(nums[1:], start=1):
        prob_new = (alpha + d*len(clusters)) / (alpha + i)
        if random.random() < prob_new:
            clusters.append([num])
            discounted_sizes.append(1 - d)
            continue

        idx = random.choices(range(len(clusters)), weights=discounted_sizes, k=1)[0]
        clusters[idx].append(num)
        discounted_sizes[idx] += 1
    return clusters


def frag(clusters: list[list[int]], d: float) -> list[list[int]]:
    fragmented = []
    for cluster in clusters:
        fragmented += crp(cluster, alpha=0, d=d)
    return fragmented


def coag(clusters: list[list[int]], alpha: float, d: float) -> list[list[int]]:
    meta_clusters = crp(range(len(clusters)), alpha=alpha/d, d=0)
    coagulated = []
    for cluster_idxs in meta_clusters:
        meta_cluster = []
        for i in cluster_idxs:
            meta_cluster += clusters[i]
        coagulated.append(meta_cluster)
    return coagulated


def generate(HP: HyperParams, seed: int = 0) -> tuple[np.ndarray, list[list[list[int]]], list[list[list[int]]]]:
    if seed != 0:
        np.random.seed(seed)
        random.seed(seed)

    # Params.
    alpha = gamma.rvs(a=HP.tau_1, scale=1/HP.tau_2)
    d = beta.rvs(a=HP.v_1, b=HP.v_2, size=HP.L-1)

    gammal = gamma.rvs(a=HP.phi_1, scale=1/HP.phi_2, size=HP.L)
    betal = np.concat([dirichlet.rvs([gammal[l]]*HP.K, size=1) for l in range(HP.L)], axis=0)

    # Cluster assignments.
    r0 = crp(range(HP.N), alpha, d=0)
    rs = [r0]
    qs = []
    for l in range(HP.L-1):
        q = frag(rs[-1], d[l])
        r = coag(q, alpha, d[l])
        qs.append(q)
        rs.append(r)

    # Cluster emissions.
    x = np.zeros(shape=(HP.N, HP.L), dtype=np.uint8)
    for l, r in enumerate(rs):
        for cluster in r:
            emission = np.random.choice(a=HP.K, p=betal[l])
            x[cluster, l] = emission
    return x, rs, qs


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--seed", default=0, type=int)
    HP = HyperParams().add_params(parser)
    args = parser.parse_args()
    HP.override_args(args)

    x, rs, qs = generate(HP, args.seed)
    print(f"r0: {rs[0]}")
    for l in range(HP.L-1):
        print(f"q{l}: {qs[l]}")
        print(f"r{l+1}: {rs[l+1]}")
    print(f"\nx:\n{x}")

