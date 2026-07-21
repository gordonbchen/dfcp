import subprocess
import json
from argparse import ArgumentParser
from functools import partial
from bayes_opt import BayesianOptimization, acquisition


def dfcp_reparam(
    seq_file: str,
    mask: float, val: float,
    tau_mu: float, tau_1: float, phi_mu: float, phi_1: float, v_mu: float, v_conc: float,
    soft: int, block_init: int
):
    tau_2 = tau_1 / tau_mu
    phi_2 = phi_1 / phi_mu

    v_1 = v_mu * v_conc
    v_2 = (1-v_mu) * v_conc

    return dfcp(
        seq_file=seq_file,
        mask=mask, val=val,
        tau_1=tau_1, tau_2=tau_2, v_1=v_1, v_2=v_2, phi_1=phi_1, phi_2=phi_2,
        soft=soft, block_init=block_init
    )


def dfcp(
    seq_file: str,
    mask: float, val: float,
    tau_1: float, tau_2: float, v_1: float, v_2: float, phi_1: float, phi_2: float,
    soft: int, block_init: int
) -> float:
    cmd = [
        "./build/dfcp", seq_file, 

        "--mask", mask, "--val", val,

        "--tau_1", str(tau_1), "--tau_2", str(tau_2),
        "--v_1", str(v_1), "--v_2", str(v_2),
        "--phi_1", str(phi_1), "--phi_2", str(phi_2),

        "--soft", str(soft),
        "--block_init", str(block_init),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        return 0.0

    out = json.loads(res.stdout)
    impute_acc = out["dfcp_impute_acc"]
    return impute_acc


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--seq_file", default="data/sims/haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt")
    p.add_argument("--mask", default="0.2")
    p.add_argument("--val", default="0.2")
    p.add_argument("--soft", type=int, default=0)
    p.add_argument("--block_init", type=int, default=0)
    p.add_argument("--init_points", type=int, default=10)
    p.add_argument("--n_iter", type=int, default=100)
    args = p.parse_args()

    f = partial(
        dfcp_reparam,
        seq_file=args.seq_file,
        mask=args.mask, val=args.val,
        soft=args.soft, block_init=args.block_init,
    )

    pbounds = {
        # beta distribution: https://www.desmos.com/calculator/bdftcfvdus
        # mu = alpha / (alpha + beta)
        # conc = alpha + beta
        "v_mu": (0, 1), "v_conc": (0.1, 200),

        # gamma distribution: https://www.desmos.com/calculator/49vz26bvex
        # mu = alpha / beta
        "tau_mu": (0.1, 100), "tau_1": (0, 100),
        "phi_mu": (0.1, 100), "phi_1": (0, 100),
    }

    optim = BayesianOptimization(
        f=f,
        acquisition_function=acquisition.UpperConfidenceBound(),
        pbounds=pbounds,
        verbose=2,
        random_state=42
    )

    subprocess.run(["./build.sh"], check=True)
    optim.maximize(init_points=args.init_points, n_iter=args.n_iter)
    print(optim.max)

    optim.save_state("optim.json")

