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

    def __eq__(self, other):
        return isinstance(other, Cluster) and (self.uuid == other.uuid)

    def __repr__(self):
        return f"E {self.emission}   {self.seqs}"

    def add_child(self, child: Cluster):
        self.children.add(child)
        child.parents.add(self)

    def remove(self, idx: int):
        self.seqs.remove(idx)
        if len(self.seqs) == 0:
            for parent in self.parents:
                parent.children.remove(self)
            for child in self.children:
                child.parents.remove(self)
            if self.emission:
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
    r = Cluster(seqs=set(range(HP.N)), emission=x_modes[0],
                seg_group=segregated_rs[0], cluster_group=rs[0])
    qs = [set() for _ in range(HP.L-1)]
    for l in range(HP.L-1):
        q = Cluster(seqs=set(range(HP.N)), cluster_group=qs[l])
        r.add_child(q)
        r = Cluster(seqs=set(range(HP.N)), emission=x_modes[l+1],
                    seg_group=segregated_rs[l+1], cluster_group=rs[l+1])
        q.add_child(r)
    pprint(rs)
    pprint(qs)

    print("segregated_rs:")
    pprint(segregated_rs)

    r_assignments = [[] for i in range(HP.N)]
    for r in rs:
        for a in r:
            for i in a.seqs:
                r_assignments[i].append(a)
    q_assignments = [[] for i in range(HP.N)]
    for q in qs:
        for b in q:
            for i in b.seqs:
                q_assignments[i].append(b)
    print("r_assignments:")
    pprint(r_assignments)
    print("q_assignments:")
    pprint(q_assignments)

    for i in range(HP.N):
        for a in r_assignments[i]:
            a.remove(i)
        for b in q_assignments[i]:
            b.remove(i)
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

            ma["new"] = (gammal[l].expect(lambda g: np.log(g + nkl) - np.log(HP.K * g + nRl)), None)
            for a in rs[l]:
                ma[a] = (0 if x[i, l] == a.emission else float("-infinity"), None)
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
            print(mb)

            ma["new"] = (ma["new"][0] + mb["new"][0], "new")
            for a in rs[l]:
                next_mb = {}
                next_mb["new"] = np.log(len(a.children)) + d[l].expect(np.log) + mb["new"][0]
                for b in a.children:
                    next_mb[b] = d[l].expect(lambda x: np.log(len(b.seqs) - x)) + mb[b][0]
                next_b = max(next_mb, key=next_mb.get)
                ma[a] = (ma[a][0] - np.log(len(a.seqs)) + next_mb[next_b], next_b)
            print(ma)

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
        print(path)

        for j in range(len(path)):
            l = j // 2
            print(j, l)
            if path[j] != "new":
                path[j].seqs.add(i)
                continue

            if j % 2 == 0:
                new_cluster = Cluster(seqs=set([i]), emission=x[i, l],
                            seg_group=segregated_rs[l], cluster_group=rs[l])
            else:
                print(qs)
                new_cluster = Cluster(seqs=set([i]), cluster_group=qs[l])
            if j - 1 >= 0:
                path[j-1].add_child(new_cluster)
            if j + 1 < len(path) and path[j+1] != "new":
                new_cluster.add_child(path[j+1])
            path[j] = new_cluster
        print(path)
        break




if __name__ == "__main__":
    # TODO: different HP for generate and me.
    HP = HyperParams().cli()
    x, _, _ = generate(HP)
    print(f"x:\n{x}")

    me(x, HP)
