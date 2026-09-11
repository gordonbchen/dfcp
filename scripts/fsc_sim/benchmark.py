#!/usr/bin/env python3

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path


FIELDS = [
    "run", "pbwt_match_len", "stage", "tau_1", "tau_2", "v_1", "v_2", "wall_s", "t_init_ms",
    "steps", "t_max_ms", "mean_nR", "mean_adj_iou", "mean_clusters", "mean_excess_parsimony",
    "clade_iou", "mean_tract_bp",
]


def evaluate(root: Path, output: Path, name: str, log_file: Path) -> dict[str, object]:
    ref = root / "data/fsc/prepared/ref.bin"
    positions = root / "data/fsc/prepared/variant_pos.txt"
    trees = root / "data/fsc/ex_0_pop_1/ex_0_pop_1_1_true_trees.trees"
    r_assign = output / f"{name}.r_assign.bin"
    eval_file = output / f"{name}.eval.json"
    clade_times = output / f"{name}.clade_times.tsv"
    cluster_tracts = output / f"{name}.cluster_tracts.tsv"
    command = [
        str(root / "build/eval_clusters"), str(ref), str(r_assign), str(positions), str(trees),
        "--clade_times", str(clade_times), "--cluster_tracts", str(cluster_tracts),
    ]
    with eval_file.open("w") as stdout, log_file.open("a") as stderr:
        subprocess.run(command, cwd=root, stdout=stdout, stderr=stderr, check=True)
    return json.loads(eval_file.read_text())


def run_one(root: Path, output: Path, name: str, match_len: int, stage: str,
            tau_1: float = 1.0, tau_2: float = 1.0,
            v_1: float = 1.0, v_2: float = 1.0) -> dict[str, object]:
    ref = root / "data/fsc/prepared/ref.bin"
    r_assign = output / f"{name}.r_assign.bin"
    result_file = output / f"{name}.json"
    log_file = output / f"{name}.log"

    command = [
        str(root / "build/impute"), str(ref), "--init", "pbwt", "--pbwt_match_len", str(match_len),
        "--tau_1", str(tau_1), "--tau_2", str(tau_2), "--v_1", str(v_1), "--v_2", str(v_2),
        "--output_r_assign", str(r_assign),
    ]
    if stage == "init":
        command += ["--init_only", "1"]
    elif stage == "step1":
        command += ["--max_train_steps", "1"]
    elif stage == "step3":
        command += ["--max_train_steps", "3"]

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    start = time.perf_counter()
    with result_file.open("w") as stdout, log_file.open("w") as stderr:
        subprocess.run(command, cwd=root, env=env, stdout=stdout, stderr=stderr, check=True)
    wall_s = time.perf_counter() - start
    fit = json.loads(result_file.read_text())

    evaluation = evaluate(root, output, name, log_file)

    train_log = fit.get("train_log", [])
    return {
        "run": name,
        "pbwt_match_len": match_len,
        "stage": stage,
        "tau_1": tau_1,
        "tau_2": tau_2,
        "v_1": v_1,
        "v_2": v_2,
        "wall_s": f"{wall_s:.3f}",
        "t_init_ms": fit["t_init"],
        "steps": len(train_log),
        "t_max_ms": sum(step["t_max"] for step in train_log),
        "mean_nR": train_log[-1]["mean_nR"] if train_log else evaluation["mean_clusters"],
        "mean_adj_iou": evaluation["mean_adj_iou"],
        "mean_clusters": evaluation["mean_clusters"],
        "mean_excess_parsimony": evaluation["mean_excess_parsimony"],
        "clade_iou": evaluation["clade_iou"],
        "mean_tract_bp": evaluation["mean_cluster_tract_span_bp"],
    }


def run_beagle4(root: Path, output: Path, javac: str, java: str) -> dict[str, object]:
    name = "beagle4"
    gen = root / "data/fsc/ex_0_pop_1/ex_0_pop_1_1_1.gen"
    log_file = output / f"{name}.log"
    command = [
        sys.executable, str(root / "scripts/fsc_sim/beagle4.py"), str(gen), str(output),
        "--javac", javac, "--java", java,
    ]
    start = time.perf_counter()
    with log_file.open("w") as stream:
        subprocess.run(command, cwd=root, stdout=stream, stderr=subprocess.STDOUT, check=True)
    evaluation = evaluate(root, output, name, log_file)
    return {
        "run": name,
        "pbwt_match_len": "",
        "stage": "dag",
        "tau_1": "",
        "tau_2": "",
        "v_1": "",
        "v_2": "",
        "wall_s": f"{time.perf_counter() - start:.3f}",
        "t_init_ms": "",
        "steps": "",
        "t_max_ms": "",
        "mean_nR": evaluation["mean_clusters"],
        "mean_adj_iou": evaluation["mean_adj_iou"],
        "mean_clusters": evaluation["mean_clusters"],
        "mean_excess_parsimony": evaluation["mean_excess_parsimony"],
        "clade_iou": evaluation["clade_iou"],
        "mean_tract_bp": evaluation["mean_cluster_tract_span_bp"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark DFCP cluster recovery on the fastsimcoal fixture."
    )
    parser.add_argument("--output", type=Path, default=Path("output/fsc_clusters"))
    parser.add_argument("--javac", default="javac")
    parser.add_argument("--java", default="java")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    output = (root / args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)

    default = (1.0, 1.0, 1.0, 1.0)
    low_d = (1.0, 1.0, 100.0, 999900.0)
    high_alpha = (1000.0, 1.0, 1.0, 1.0)
    low_d_high_alpha = (1000.0, 1.0, 100.0, 999900.0)
    runs = [
        ("pbwt50_init", 50, "init", *default),
        ("pbwt100_init", 100, "init", *default),
        ("pbwt200_init", 200, "init", *default),
        ("pbwt200_step1", 200, "step1", *default),
        ("pbwt200_step3", 200, "step3", *default),
        ("pbwt200_converged", 200, "converged", *default),
        ("pbwt200_low_d_step3", 200, "step3", *low_d),
        ("pbwt200_high_alpha_step3", 200, "step3", *high_alpha),
        ("pbwt200_low_d_high_alpha_step3", 200, "step3", *low_d_high_alpha),
    ]

    summary = output / "summary.tsv"
    with summary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t")
        writer.writeheader()
        for name, match_len, stage, tau_1, tau_2, v_1, v_2 in runs:
            print(f"running {name}", flush=True)
            row = run_one(root, output, name, match_len, stage, tau_1, tau_2, v_1, v_2)
            writer.writerow(row)
            stream.flush()
        print("running beagle4", flush=True)
        writer.writerow(run_beagle4(root, output, args.javac, args.java))
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
