import random
from dataclasses import dataclass

import numpy as np
from scipy.stats import gamma, beta, dirichlet

from cli_params import CLIParams


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


def generate(HP: HyperParams) -> tuple[np.ndarray, list[list[int]], list[list[int]]]:
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
    HP = HyperParams().cli()

    x, rs, qs = generate(HP)
    print(f"r0: {rs[0]}")
    for l in range(HP.L-1):
        print(f"q{l}: {qs[l]}")
        print(f"r{l+1}: {rs[l+1]}")
    print(f"\nx:\n{x}")

