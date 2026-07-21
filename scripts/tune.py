import subprocess
import json
from argparse import ArgumentParser
from functools import partial
from bayes_opt import BayesianOptimization, acquisition


def get_dfcp_cmd(
    tau_1: float, tau_2: float, v_1: float, v_2: float, phi_1: float, phi_2: float,
    seq_file: str,
    mask: float, val: float,
    soft: int, block_init: int,
    tree_fname: str | None, variant_pos_fname: str | None,
    variant_start_pos: int | None, clade_beta: float | None,
) -> list[str]:
    cmd = [
        "./build/dfcp", seq_file, 

        "--mask", str(mask), "--val", str(val),

        "--tau_1", str(tau_1), "--tau_2", str(tau_2),
        "--v_1", str(v_1), "--v_2", str(v_2),
        "--phi_1", str(phi_1), "--phi_2", str(phi_2),

        "--soft", str(soft), "--block_init", str(block_init)
    ]
    if not tree_fname is None:
        assert not None in (variant_pos_fname, variant_start_pos, clade_beta)
        cmd += [
            "--tree", tree_fname,
            "--variant_pos_fname", variant_pos_fname,
            "--variant_start_pos", str(variant_start_pos),
            "--clade_beta", str(clade_beta),
        ]
    return cmd


def dfcp(
    tau_1: float, tau_2: float, v_1: float, v_2: float, phi_1: float, phi_2: float,
    opt_name: str,
    seq_file: str,
    mask: float, val: float,
    soft: int, block_init: int,
    tree_fname: str | None, variant_pos_fname: str | None,
    variant_start_pos: int | None, clade_beta: float | None,
) -> float:
    cmd = get_dfcp_cmd(
        tau_1=tau_1, tau_2=tau_2, v_1=v_1, v_2=v_2, phi_1=phi_1, phi_2=phi_2,
        seq_file=seq_file,
        mask=mask, val=val,
        soft=soft, block_init=block_init,
        tree_fname=tree_fname, variant_pos_fname=variant_pos_fname,
        variant_start_pos=variant_start_pos, clade_beta=clade_beta,
    )
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
        soft=args.soft, block_init=args.block_init,
        tree_fname=args.tree_fname, variant_pos_fname=args.variant_pos_fname,
        variant_start_pos=args.variant_start_pos, clade_beta=args.clade_beta,
    )
    f = partial(
        dfcp,
        opt_name=args.opt_name,
        **eval_args
    )

    pbounds = dict(
        tau_1=(0.01, 10), tau_2=(0.01, 10),
        v_1=(0.01, 10), v_2=(0.01, 10),
        phi_1=(0.01, 10), phi_2=(0.01, 10),
    )

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

