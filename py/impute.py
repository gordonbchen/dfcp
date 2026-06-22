from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from generate import generate, HyperParams, seed_rngs
from me import me


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--seed", default=0, type=int)
    parser.add_argument("--seq_file", default=None)
    HP = HyperParams().add_params(parser)
    args = parser.parse_args()
    seed_rngs(args.seed)
    HP.override_args(args)

    if args.seq_file is None:
        x, _, _ = generate(HP)
    else:
        with open(args.seq_file, "r") as f:
            x = np.array([list(map(int, line)) for line in f.read().split()], dtype=np.int8)
        HP.K = int(np.max(x)) + 1
        HP.N, HP.L = x.shape
    print(HP)
    print(f"x:\n{x}")
    x_true = x.copy()

    mask_fracs = [0.001, 0.01, 0.1, 0.2, 0.5]
    accs = []
    baseline_accs = []
    for mask_frac in mask_fracs:
        print(f"\nmask_frac={mask_frac}")
        x = x_true.copy()
        x_mask = np.random.random(x.shape) < mask_frac
        x[x_mask] = -1

        total_masked = x_mask.sum()
        baseline_correct = 0
        for l in range(HP.L):
            m = stats.mode(x[:, l][x[:, l] >= 0]).mode
            baseline_correct += (x_true[:, l][x_mask[:, l]] == m).sum()
        baseline_acc = baseline_correct / total_masked
        baseline_accs.append(baseline_acc)
        print(f"baseline accuracy: {baseline_correct} / {total_masked} = {baseline_acc}")

        r_assignments = me(x, HP)

        correct = 0
        for i in range(HP.N):
            for l in range(HP.L):
                if x_mask[i, l] and r_assignments[i][l].emission == x_true[i, l]:
                    correct += 1
        acc = correct / total_masked
        accs.append(acc)
        print(f"imputation accuracy: {correct} / {total_masked} = {acc}")

    plt.plot(mask_fracs, accs, "-bo", label="dfcp")
    plt.plot(mask_fracs, baseline_accs, "-ro", label="loc mode baseline")
    plt.title("Imputation Accuracy")
    plt.xlabel("mask frac")
    plt.ylabel("imputation acc")
    plt.legend(loc=0)
    plt.savefig("impute_acc.png")
    plt.close()

