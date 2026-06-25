import subprocess
import re
import matplotlib.pyplot as plt
from pathlib import Path


xs, ys, parsimonys = [], [], []

r_file = re.compile(r"haps_SIMOUT_1.txt.gz_(\d+\.\d+)_(\d+\.\d+)")
r_parsimony = re.compile(r"mean_excess_parsimony=(\d+\.\d+)")

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

    m = r_parsimony.search(res.stdout)
    parsimony = float(m.group(1))

    print(x, y, parsimony)

    xs.append(x)
    ys.append(y)
    parsimonys.append(parsimony)


for y0 in sorted(set(ys)):
    xs_y, parsimonys_y = [], []
    for i, y in enumerate(ys):
        if y == y0:
            xs_y.append(xs[i])
            parsimonys_y.append(parsimonys[i])

    idxs = sorted(range(len(xs_y)), key=lambda i: xs_y[i])
    xs_y = [xs_y[i] for i in idxs]
    parsimonys_y = [parsimonys_y[i] for i in idxs]

    plt.plot(xs_y, parsimonys_y, "-o", label=f"y={y0}")
plt.ylabel("mean excess parsimony")
plt.xlabel("x")
plt.xscale("symlog", linthresh=0.001)
plt.legend(loc=0)
plt.savefig("parsimony.png")

