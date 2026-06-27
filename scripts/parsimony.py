import subprocess
import re
import json
import matplotlib.pyplot as plt
from pathlib import Path
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument("--out_file", default="parsimony.png", required=False)
args = parser.parse_args()

xs, ys, parsimonys, emission_parsimonys = [], [], [], []

r_file = re.compile(r"haps_SIMOUT_1.txt.gz_(\d+\.\d+)_(\d+\.\d+)")

# Run inference.
seq_dir = Path("data/examples/simulated/SIM1_LEN500_NHAPS100")
for seq_file in seq_dir.glob("haps_*.txt"):
    m = r_file.match(seq_file.name)
    x, y = (0.0, 0.0) if m is None else map(float, m.groups())

    res = subprocess.run([
        "./build/dfcp", seq_file,
        "--tree", "data/examples/simulated/SIM1_LEN500_NHAPS100/ex_0_pop_1_1_true_trees.trees_1_14572_500.trees",
        "--variant_pos", "data/examples/simulated/SIM1_LEN500_NHAPS100/variant_pos.txt",
        "--variant_start_pos", "14572"
    ], capture_output=True, text=True)

    out = json.loads(res.stdout)
    parsimony = float(out["mean_excess_parsimony"])
    emission_parsimony = float(out["mean_emission_excess_parsimony"])
    print(x, y, parsimony, emission_parsimony)

    xs.append(x)
    ys.append(y)
    parsimonys.append(parsimony)
    emission_parsimonys.append(emission_parsimony)


# Plot.
for ps, name in zip((parsimonys, emission_parsimonys), ("dfcp", "emission")):
    for y0 in sorted(set(ys)):
        xs_y, parsimonys_y = [], []
        for i, y in enumerate(ys):
            if y == y0:
                xs_y.append(xs[i])
                parsimonys_y.append(ps[i])

        idxs = sorted(range(len(xs_y)), key=lambda i: xs_y[i])
        xs_y = [xs_y[i] for i in idxs]
        parsimonys_y = [parsimonys_y[i] for i in idxs]

        shape = "x" if name == "dfcp" else "o"
        plt.plot(xs_y, parsimonys_y, f"-{shape}", label=f"{name}: y={y0}")

plt.ylabel("mean excess parsimony")
plt.xlabel("x")
plt.xscale("symlog", linthresh=0.001)
plt.legend(loc=0)
plt.savefig(args.out_file)

