from argparse import ArgumentParser
import numpy as np


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--hap_file", default="pbwt_test.txt")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    with open(args.hap_file, "r") as f:
        text = f.read()
    X = np.array([[int(x) for x in l] for l in text.split()], dtype=np.int8)
    N, L = X.shape
    if args.verbose: print(f"X {X.shape}\n{X}")

    # Build prefix and divergence arrays.
    a = np.zeros((N,L+1), dtype=np.int32)
    a[:, 0] = np.arange(N)

    d = np.zeros((N,L+1), dtype=np.int32)

    for l in range(L):
        zeros, ones = [], []

        zero_d, one_d = [], []
        p, q = l+1, l+1

        for i in range(N):
            p = max(p, d[i,l])
            q = max(q, d[i,l])

            if X[a[i,l],l] == 0:
                zeros.append(a[i,l])
                zero_d.append(p)
                p = 0
            else:
                ones.append(a[i,l])
                one_d.append(q)
                q = 0

        a[:,l+1] = zeros + ones
        d[:N,l+1] = zero_d + one_d
    if args.verbose: print(f"a\n{a}\nd\n{d}")

    # PBWTed X.
    Y = np.zeros((N,L), np.int32)
    for l in range(L):
        for i in range(N):
            Y[i,l] = X[a[i,l],l]
    if args.verbose: print(f"Y\n{Y}")

    # Locally maximal matches of length at least K.
    # A prefix at position l starting at d[i,l] has length l-d[i,l].
    # For a match of at least K: l-d[i,l] >= K, so d[i,l] <= l-K.
    K = 2
    print(f"Finding locally maximal matches of length at least {K}.")
    local_matches = []
    for l in range(K, L):
        zeros, ones = [], []

        for i in range(N):
            if d[i,l] > l-K:
                for j in zeros:
                    for k in ones:
                        local_matches.append((l, j, k))
                zeros.clear()
                ones.clear()

            if Y[i,l] == 0:
                zeros.append(int(a[i,l]))
            else:
                ones.append(int(a[i,l]))

        for j in zeros:
            for k in ones:
                local_matches.append((l, j, k))

    match_group = []  # Final location matches all combinations, no zeros x ones.
    for i in range(N):
        if d[i,L] > L-K:
            for j in range(len(match_group)):
                for k in range(j+1, len(match_group)):
                    local_matches.append((L, match_group[j], match_group[k]))
            match_group.clear()
        match_group.append(int(a[i,L]))
    for j in range(len(match_group)):
        for k in range(j+1, len(match_group)):
            local_matches.append((L, match_group[j], match_group[k]))

    if args.verbose: print(f"local matches: {local_matches}")
    print(f"# local matches: {len(local_matches)}")

    # Set maximal matches.
    set_max_matches = []
    for l in range(1, L+1):
        for i in range(N):
            lb, ub = i, i
            loc_match = False

            # Extend back.
            if (i == N-1) or ((i != 0) and (d[i,l] <= d[i+1,l])):
                while (lb > 0) and (d[lb,l] <= d[i,l]):
                    lb -= 1
                    # Not set maximal if can be right extended.
                    if (l != L) and (Y[i,l] == Y[lb,l]):
                        loc_match = True
                        break
                if loc_match: continue

            # Extend forwards.
            if (i != N-1) and (d[i+1,l] <= d[i,l]):
                while (ub < N-1) and (d[ub+1,l] <= d[i+1,l]):
                    ub += 1
                    if (l != L) and (Y[i,l] == Y[ub,l]):
                        loc_match = True
                        break
                if loc_match: continue

            for j in range(lb, i):
                if d[i,l] < l:
                    set_max_matches.append((int(d[i,l]), l, int(a[i,l]), int(a[j,l])))
            for j in range(i+1, ub+1):
                if d[i+1,l] < l:
                    set_max_matches.append((int(d[i+1,l]), l, int(a[i,l]), int(a[j,l])))
    if args.verbose: print(f"set maximal matches: {set_max_matches}")
    print(f"# set maximal matches: {len(set_max_matches)}")

    # Finding maximal matches for a new haplotype.
    z = np.array([1, 1, 1, 0, 0], np.int32)
    if args.verbose: print(f"z\n{z}")

    u = np.zeros((N+1,L), np.int32)
    v = np.zeros((N+1,L), np.int32)

    for l in range(L):
        for i in range(N):
            u[i+1,l] = u[i,l]
            v[i+1,l] = v[i,l]

            if Y[i,l] == 0:
                u[i+1,l] += 1
            else:
                v[i+1,l] += 1
    if args.verbose: print(f"u\n{u}\nv\n{v}")

    matches = []
    e, f, g = 0, 0, N
    for l in range(L):
        f_new = u[f,l] if z[l] == 0 else u[N,l] + v[f,l]
        g_new = u[g,l] if z[l] == 0 else u[N,l] + v[g,l]

        if f_new < g_new:
            f, g = f_new, g_new
            continue

        if e < l:
            for i in range(f, g):
                matches.append((e, l, int(a[i,l])))

        f, g = f_new, g_new

        e = l+1 if (f == N) or (f == 0) else int(d[f,l+1]) - 1
        if (f == N) or ((f > 0) and (z[e] == 0)):
            # Extend back.
            f -= 1
            while (e > 0) and (z[e-1] == X[a[f,l+1],e-1]): e -= 1
            while (f > 0) and (d[f,l+1] <= e): f -= 1
        else:
            # Extend forwards.
            g += 1
            while (e > 0) and (z[e-1] == X[a[f,l+1],e-1]): e -= 1
            while (g < N) and (d[g,l+1] <= e): g += 1

    if e < L:
        for i in range(f, g):
            matches.append((e, L, int(a[i,L])))

    print(f"# set maximal matches to z: {len(matches)}")
    if args.verbose: print(f"set maximal matches to z: {matches}")

