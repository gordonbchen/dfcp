"""Shared subprocess interface to the C++ DFCP executable."""

import json
import subprocess
from argparse import ArgumentParser
from pathlib import Path


def build_dfcp() -> None:
    subprocess.run(["./build.sh"], check=True)


def get_dfcp_cmd(seq_file: str | Path, **options: object) -> list[str]:
    command = ["./build/dfcp", str(seq_file)]
    for name, value in options.items():
        if value is not None:
            command.extend((f"--{name}", str(value)))
    return command


def run_dfcp(seq_file: str | Path, retries: int = 0, **options: object) -> dict:
    command = get_dfcp_cmd(seq_file, **options)
    last_error = ""
    for _ in range(retries + 1):
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            try:
                return json.loads(result.stdout)
            except json.JSONDecodeError as error:
                last_error = f"DFCP emitted invalid JSON: {error}"
        else:
            last_error = result.stderr.strip()[-1000:]
    raise RuntimeError(last_error or "DFCP failed without an error message")


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

    p.add_argument("--mode", choices=("hard", "noisy", "soft"), default=None)
    p.add_argument("--lambda_1", type=float, default=None)
    p.add_argument("--lambda_2", type=float, default=None)

    p.add_argument("--init", choices=("viterbi", "block", "pbwt"), default=None)
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

        mode=args.mode,
        lambda_1=args.lambda_1, lambda_2=args.lambda_2,

        init=args.init, pbwt_match_len=args.pbwt_match_len,
        init_only=args.init_only,

        tau_1=args.tau_1, tau_2=args.tau_2,
        v_1=args.v_1, v_2=args.v_2,
        phi_1=args.phi_1, phi_2=args.phi_2,
    )

    build_dfcp()
    res = run_dfcp(args.seq_file, **dfcp_kwargs)
    print(res)
