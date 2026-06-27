import subprocess
import random
import math
import json
from argparse import ArgumentParser


p = ArgumentParser()
p.add_argument("--seq_file", default="data/sims/haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt")
p.add_argument("--mask", default="0.2")
p.add_argument("--val", default="0.2")
p.add_argument("--trials", type=int, default=100)
p.add_argument("--seed", type=int, default=0)
args = p.parse_args()


rng = random.Random(args.seed)
logu = lambda lo, hi: math.exp(rng.uniform(math.log(lo), math.log(hi)))

def gamma_hp(mean_lo=0.01, mean_hi=100, shape_lo=0.05, shape_hi=50):
    mean, shape = logu(mean_lo, mean_hi), logu(shape_lo, shape_hi)
    return shape, shape / mean

def beta_hp(mean_lo=0.01, mean_hi=0.99, conc_lo=0.1, conc_hi=1000):
    m, k = rng.uniform(mean_lo, mean_hi), logu(conc_lo, conc_hi)
    return m * k, (1 - m) * k


subprocess.run(["./build.sh"], check=True)
best = None

for i in range(args.trials):
    tau_1, tau_2 = gamma_hp()
    phi_1, phi_2 = gamma_hp()
    v_1, v_2 = beta_hp()

    cmd = ["./build/dfcp", args.seq_file, "--mask", args.mask, "--val", args.val,
           "--tau_1", str(tau_1), "--tau_2", str(tau_2),
           "--v_1", str(v_1), "--v_2", str(v_2),
           "--phi_1", str(phi_1), "--phi_2", str(phi_2)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        print(f"{i}: FAILED", res.stderr.strip() or res.stdout[-500:], "\n")
        continue

    out = json.loads(res.stdout)
    acc = out["dfcp_impute_acc"]
    elbo = max(d["elbo"] for d in out["train_log"])

    cur = dict(acc=acc, elbo=elbo, tau_1=tau_1, tau_2=tau_2, v_1=v_1, v_2=v_2, phi_1=phi_1, phi_2=phi_2)
    if best is None or (acc, elbo) > (best["acc"], best["elbo"]):
        best = cur

    print(f"{i}: acc={acc:.4f} elbo={elbo:.4f} "
          f"tau=({tau_1:.4g},{tau_2:.4g}) v=({v_1:.4g},{v_2:.4g}) phi=({phi_1:.4g},{phi_2:.4g}) "
          f"best_acc={best['acc']:.4f}")

print("BEST:", best)

