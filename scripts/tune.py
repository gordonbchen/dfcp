import subprocess
import json
from argparse import ArgumentParser
from functools import partial
from bayes_opt import BayesianOptimization, acquisition


def get_dfcp_cmd(
    seq_file: str,
    mask: float, val: float,
    tau_1: float, tau_2: float, v_1: float, v_2: float, phi_1: float, phi_2: float,
    soft: int, block_init: int
) -> list[str]:
    return [
        "./build/dfcp", seq_file, 

        "--mask", str(mask), "--val", str(val),

        "--tau_1", str(tau_1), "--tau_2", str(tau_2),
        "--v_1", str(v_1), "--v_2", str(v_2),
        "--phi_1", str(phi_1), "--phi_2", str(phi_2),

        "--soft", str(soft),
        "--block_init", str(block_init),
    ]


def dfcp(
    seq_file: str,
    mask: float, val: float,
    tau_1: float, tau_2: float, v_1: float, v_2: float, phi_1: float, phi_2: float,
    soft: int, block_init: int
) -> float:
    cmd = get_dfcp_cmd(
        seq_file=seq_file,
        mask=mask, val=val,
        tau_1=tau_1, tau_2=tau_2,
        v_1=v_1, v_2=v_2,
        phi_1=phi_1, phi_2=phi_2,
        soft=soft, block_init=block_init
    )
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
        dfcp,
        seq_file=args.seq_file,
        mask=args.mask, val=args.val,
        soft=args.soft, block_init=args.block_init,
    )

    pbounds = {
        "tau_1": (0.01, 10), "tau_2": (0.01, 10),
        "v_1": (0.01, 10), "v_2": (0.01, 10),
        "phi_1": (0.01, 10), "phi_2": (0.01, 10),
    }

    optim = BayesianOptimization(
        f=f,
        acquisition_function=acquisition.UpperConfidenceBound(),
        pbounds=pbounds,
        verbose=2,
    )

    subprocess.run(["./build.sh"], check=True)
    optim.maximize(init_points=args.init_points, n_iter=args.n_iter)
    optim.save_state("optim.json")

    print(optim.max, "\n")
    cmd = get_dfcp_cmd(
        seq_file=args.seq_file,
        mask=args.mask, val=args.val,
        soft=args.soft, block_init=args.block_init,
        **optim.max["params"]
    )
    print(" ".join(cmd))

