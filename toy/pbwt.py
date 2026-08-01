from argparse import ArgumentParser
import numpy as np


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--hap_file", default="data/examples/simulated/SIM1_LEN500_NHAPS100/" \
                                              "haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt")
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

    d = np.zeros((N+1,L+1), dtype=np.int32)
    d[-1] = np.arange(L+1)

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

            if X[a[i,l],l] == 0:
                zeros.append(a[i,l])
            else:
                ones.append(a[i,l])

        for j in zeros:
            for k in ones:
                local_matches.append((l, j, k))

    match_group = []  # Final location matches all combinations, no zeros x ones.
    for i in range(N):
        if d[i,L] > L-K:
            for j in range(len(match_group)):
                for k in range(j+1, len(match_group)):
                    local_matches.append((L, j, k))
            match_group.clear()
        match_group.append(a[i,L])
    for j in range(len(match_group)):
        for k in range(j+1, len(match_group)):
            local_matches.append((L, j, k))

    if args.verbose: print(f"local matches: {local_matches}")
    print(f"# local matches: {len(local_matches)}")

    # Set maximal matches.
    set_max_matches = []
    for l in range(1, L+1):
        for i in range(N):
            lb, ub = i, i
            loc_match = False

            # Extend back.
            if (i != 0) and (d[i,l] <= d[i+1,l]):
                while (lb > 0) and (d[lb,l] <= d[i,l]):
                    lb -= 1
                    # Not set maximal if can be right extended.
                    if (l != L) and (X[a[i,l],l] == X[a[lb,l],l]):
                        loc_match = True
                        break
                if loc_match: continue

            # Extend forwards.
            if (i != N-1) and (d[i+1,l] <= d[i,l]):
                while (ub < N-1) and (d[ub,l] <= d[i+1,l]):
                    ub += 1
                    if (l != L) and (X[a[i,l],l] == X[a[ub,l],l]):
                        loc_match = True
                        break
                if loc_match: continue

            for j in range(lb, i):
                set_max_matches.append((d[i,l], l, a[i,l], a[j,l]))
            for j in range(i+1, ub+1):
                set_max_matches.append((d[i+1,l], l, a[i,l], a[j,l]))
    if args.verbose: print(f"set maximal matches: {set_max_matches}")
    print(f"# set maximal matches: {len(set_max_matches)}")

