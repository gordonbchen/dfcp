from __future__ import annotations
import uuid
from pprint import pprint

import numpy as np
from scipy.stats import gamma, beta, mode

from generate import generate, HyperParams


class Cluster:
    def __init__(
        self, seqs: set[int], emission: int = None,
        seg_group: set[Cluster] = None, cluster_group: set[Cluster] = None
    ):
        self.uuid = uuid.uuid4().int

        self.seqs = seqs
        self.emission = emission

        self.seg_group = seg_group
        if not self.emission is None:
            self.seg_group[emission].add(self)

        self.cluster_group = cluster_group
        self.cluster_group.add(self)

        self.children = set()
        self.parents = set()

    def __hash__(self):
        return self.uuid

    def __eq__(self, other: Cluster):
        return self.uuid == other.uuid

    def __repr__(self):
        return f"E {self.emission}   {self.seqs}"

    def add_child(self, child: Cluster):
        self.children.add(child)
        child.parents.add(self)

    def remove(self, idx: int):
        self.seqs.remove(idx)
        if len(self.seqs) == 0:
            for parent in parents:
                parent.children.remove(self)
            for child in children:
                child.parents.remove(self)
            self.seg_group[self.emission].remove(self)
            self.cluster_group.remove(self)


def me(x: np.ndarray, HP: HyperParams):
    # Init parameters.
    alpha = gamma(a=HP.tau_1, scale=1/HP.tau_2)
    d = [beta(a=HP.v_1, b=HP.v_2) for i in range(HP.L-1)]
    gammal = [gamma(a=HP.phi_1, scale=1/HP.phi_2) for i in range(HP.L)]

    # Init clusters.
    segregated_rs = [[set() for _ in range(HP.K)] for _ in range(HP.L)]
    x_modes = mode(x, axis=0)[0]
    rs = [set() for _ in range(HP.L)]
    r = Cluster(seqs=list(range(HP.N)), emission=x_modes[0],
                seg_group=segregated_rs[0], cluster_group=rs[0])
    qs = [set() for _ in range(HP.L-1)]
    for l in range(HP.L-1):
        q = Cluster(seqs=list(range(HP.N)), cluster_group=qs[l])
        r.add_child(q)
        r = Cluster(seqs=list(range(HP.N)), emission=x_modes[l+1],
                    seg_group=segregated_rs[l+1], cluster_group=rs[l+1])
        q.add_child(r)
    pprint(rs)
    pprint(qs)

    print("segregated_rs:")
    pprint(segregated_rs)

    r_assignments = [[] for i in range(HP.N)]
    for r in rs:
        for cluster in r:
            for i in cluster.seqs:
                r_assignments[i].append(cluster)
    q_assignments = [[] for i in range(HP.N)]
    for q in qs:
        for cluster in q:
            for i in cluster.seqs:
                q_assignments[i].append(cluster)
    print("r_assignments:")
    pprint(r_assignments)
    print("q_assignments:")
    pprint(q_assignments)

    for i in range(HP.N):
        for cluster in r_assignments[i]:
            cluster.remove(i)
        for cluster in q_assignments[i]:
            cluster.remove(i)
        pprint(r_assignments[i])
        pprint(q_assignments[i])
        pprint(segregated_rs)
        pprint(rs)
        pprint(qs)

        messages = []
        for l in reversed(range(HP.L)):
            # Compute likelihood term for a-message.
            ma = {}
            nkl = len(segregated_rs[l][x[i, l]])
            nRl = len(rs[l])

            ma["new"] = gammal[l].expect(lambda g: np.log(g + nkl) - np.log(HP.K * g + nRl))
            for cluster in rs[l]:
                ma[cluster] = 0 if x[i, l] == cluster.emission else float("-infinity")
            print(x)
            print(ma)
            print(l)
            if l == HP.L - 1:
                messages.append((ma,))
                continue
            print(messages)

            mb = {}
            for b in qs[l]:
                assert len(b.children) == 1
                mb[b] = messages[-1][0][next(iter(b.children))]

            nQl = len(qs[l])
            mu_y = alpha.mean() + nQl * d[l].mean()
            sigma2_y = alpha.var() + nQl*nQl * d[l].var()
            elogy = np.log(mu_y) - 0.5 * sigma2_y / (mu_y*mu_y)
            a_messages = [] # TODO:
            mb["new"] = -elogy + max(a_messages)

            print(mb)


            messages.append((ma, mb))
            break
        break




if __name__ == "__main__":
    # TODO: different HP for generate and me.
    HP = HyperParams().cli()
    x, _, _ = generate(HP)
    print(f"x:\n{x}")

    me(x, HP)
