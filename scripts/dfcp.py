import subprocess
import json
from argparse import ArgumentParser


def build_dfcp():
    subprocess.run(["./build.sh"], check=True)


def get_dfcp_cmd(seq_file: str, **kwargs) -> list[str]:
    cmd = ["./build/dfcp", seq_file]
    for k, v in kwargs.items():
        if not v is None:
            cmd.append(f"--{k}")
            cmd.append(str(v))
    return cmd


def run_dfcp(seq_file: str, **kwargs) -> dict:
    cmd = get_dfcp_cmd(seq_file, **kwargs)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        raise RuntimeError(res.stderr)
    out = json.loads(res.stdout)
    return out


def get_dfcp_parser() -> ArgumentParser:
    p = ArgumentParser()
    p.add_argument(
        "--seq_file",
        type=str,
        default="data/examples/simulated/SIM1_LEN500_NHAPS100/haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt"
    )

    p.add_argument("--val", type=float, default=None)
    p.add_argument("--mask", type=float, default=None)

    p.add_argument("--tree", type=str, default=None)
    p.add_argument("--variant_pos_fname", type=str, default=None)
    p.add_argument("--variant_start_pos", type=int, default=None)
    p.add_argument("--clade_beta", type=float, default=None)
    p.add_argument("--tree_vis", type=str, default=None)

    p.add_argument("--noisy", type=int, default=None)
    p.add_argument("--lambda_1", type=float, default=None)
    p.add_argument("--lambda_2", type=float, default=None)

    p.add_argument("--soft", type=int, default=None)

    p.add_argument("--block_init", type=int, default=None)
    p.add_argument("--pbwt_init", type=int, default=None)
    p.add_argument("--pbwt_match_len", type=int, default=None)
    p.add_argument("--init_only", type=int, default=None)

    p.add_argument("--tau_1", type=float, default=None)
    p.add_argument("--tau_2", type=float, default=None)
    p.add_argument("--v_1", type=float, default=None)
    p.add_argument("--v_2", type=float, default=None)
    p.add_argument("--phi_1", type=float, default=None)
    p.add_argument("--phi_2", type=float, default=None)
    return p


if __name__ == "__main__":
    p = get_dfcp_parser()
    args = p.parse_args()

    dfcp_kwargs = dict(
        val=args.val, mask=args.mask,

        tree=args.tree, variant_pos_fname=args.variant_pos_fname,
        variant_start_pos=args.variant_start_pos,
        tree_vis=args.tree_vis,
        clade_beta=args.clade_beta,

        noisy=args.noisy,
        lambda_1=args.lambda_1, lambda_2=args.lambda_2,

        soft=args.soft,

        block_init=args.block_init,
        pbwt_init=args.pbwt_init, pbwt_match_len=args.pbwt_match_len,
        init_only=args.init_only,

        tau_1=args.tau_1, tau_2=args.tau_2,
        v_1=args.v_1, v_2=args.v_2,
        phi_1=args.phi_1, phi_2=args.phi_2,
    )

    build_dfcp()
    res = run_dfcp(args.seq_file, **dfcp_kwargs)
    print(res)

