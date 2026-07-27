import subprocess
import json
import math
from argparse import ArgumentParser
from functools import partial
from bayes_opt import BayesianOptimization, acquisition


def get_dfcp_cmd(
    seq_file: str,
    noisy: int, soft: int, block_init: int,
    mask: float, val: float,
    tree_fname: str, variant_pos_fname: str, variant_start_pos: int, clade_beta: float,
    log_tau_1: float, log_tau_2: float, log_v_1: float, log_v_2: float, log_phi_1: float, log_phi_2: float,
    log_lambda_1: float | None = None, log_lambda_2: float | None = None,
) -> list[str]:
    tau_1, tau_2, v_1, v_2, phi_1, phi_2 = (
        math.exp(x) for x in (log_tau_1, log_tau_2, log_v_1, log_v_2, log_phi_1, log_phi_2)
    )
    cmd = [
        "./build/dfcp", seq_file, 
        "--noisy", str(noisy), "--soft", str(soft), "--block_init", str(block_init),
        "--mask", str(mask), "--val", str(val),
    ]
    if not tree_fname is None:
        assert not None in (variant_pos_fname, variant_start_pos, clade_beta)
        cmd += [
            "--tree", tree_fname,
            "--variant_pos_fname", variant_pos_fname,
            "--variant_start_pos", str(variant_start_pos),
            "--clade_beta", str(clade_beta),
        ]
    cmd += [
        "--tau_1", str(tau_1), "--tau_2", str(tau_2),
        "--v_1", str(v_1), "--v_2", str(v_2),
        "--phi_1", str(phi_1), "--phi_2", str(phi_2),
    ]
    if noisy == 1:
        assert not (log_lambda_1 is None or log_lambda_2 is None)
        lambda_1, lambda_2 = (math.exp(x) for x in (log_lambda_1, log_lambda_2))
        cmd += ["--lambda_1", str(lambda_1), "--lambda_2", str(lambda_2)]
    return cmd


def dfcp(opt_name: str, **kwargs) -> float:
    cmd = get_dfcp_cmd(**kwargs)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        return 0.0
    out = json.loads(res.stdout)
    return out[opt_name]


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "--seq_file",
        type=str,
        default="data/examples/simulated/SIM1_LEN500_NHAPS100/haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt"
    )

    p.add_argument("--mask", type=float, default=0.2)
    p.add_argument("--val", type=float, default="0.2")

    p.add_argument("--noisy", type=int, default=0)
    p.add_argument("--soft", type=int, default=0)
    p.add_argument("--block_init", type=int, default=0)

    p.add_argument("--tree_fname", type=str, default=None)
    p.add_argument("--variant_pos_fname", type=str, default=None)
    p.add_argument("--variant_start_pos", type=int, default=None)
    p.add_argument("--clade_beta", type=float, default=None)

    p.add_argument("--init_points", type=int, default=10)
    p.add_argument("--n_iter", type=int, default=100)
    p.add_argument("--opt_name", type=str, default="dfcp_impute_acc")
    args = p.parse_args()

    eval_args = dict(
        seq_file=args.seq_file,
        mask=args.mask, val=args.val,
        noisy=args.noisy, soft=args.soft, block_init=args.block_init,
        tree_fname=args.tree_fname, variant_pos_fname=args.variant_pos_fname,
        variant_start_pos=args.variant_start_pos, clade_beta=args.clade_beta,
    )
    f = partial(
        dfcp,
        opt_name=args.opt_name,
        **eval_args
    )

    pbounds = dict(
        log_tau_1=(-4, 3), log_tau_2=(-4, 3),
        log_v_1=(-4, 3), log_v_2=(-4, 3),
        log_phi_1=(-4, 3), log_phi_2=(-4, 3),
    )
    if args.noisy == 1:
        pbounds["log_lambda_1"] = (-4, 3)
        pbounds["log_lambda_2"] = (-4, 3)

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
    cmd = "./build.sh && " + " ".join(get_dfcp_cmd(**eval_args, **optim.max["params"]))
    if not args.tree_fname is None:
        cmd += " --tree_vis output/tree.dot > /dev/null && dot -Tsvg output/tree.dot -o output/tree.svg"
    print(cmd)

