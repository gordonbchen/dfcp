#!/usr/bin/env python3

"""
Split 1000 Genomes samples into a reference panel and target individuals.
Choose targets by population. Create reference and target VCFs at reference-polymorphic sites.
"""

import argparse
import csv
import random
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def record_count(path: Path) -> int:
    res = subprocess.run(["bcftools", "index", "--nrecords", str(path)],
                         check=True, text=True, capture_output=True)
    return int(res.stdout.strip())

def parse_args() -> argparse.Namespace:
    data_dir = Path("data/1000g_phase3_v5b")
    default_vcf = data_dir / "biallelic_snps" / "chr20.phase3_v5b.biallelic_snvs.vcf.gz"
    default_panel = data_dir / "integrated_call_samples_v3.20130502.ALL.panel"

    parser = argparse.ArgumentParser(
        description=("Choose targets by population. Create reference and "
                     "target VCFs at reference-polymorphic sites.")
    )
    parser.add_argument("--vcf", type=Path, default=default_vcf)
    parser.add_argument("--panel", type=Path, default=default_panel)
    parser.add_argument("--output-dir", type=Path, default=data_dir / "ref_target")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--targets-per-pop", type=int, default=2)
    parser.add_argument("--threads", type=int, default=1)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Read panel file.
    print("Reading panel file.")
    with args.panel.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    pops = defaultdict(list)
    for row in rows:
        sample = row["sample"]
        pop = row["pop"]
        pops[pop].append(sample)

    # Check # targets per pop is reachable.
    for pop, samples in pops.items():
        if len(samples) < args.targets_per_pop:
            raise RuntimeError(f"{pop} has too few samples for desired # targets.")

    # Check vcf and panel samples are the same.
    res = subprocess.run(["bcftools", "query", "-l", str(args.vcf)],
                         check=True, text=True, capture_output=True)
    vcf_sample_list = res.stdout.splitlines()
    vcf_samples = set(vcf_sample_list)
    panel_samples = {row["sample"] for row in rows}
    if len(vcf_sample_list) != len(vcf_samples):
        raise RuntimeError("VCF contains duplicate sample names.")
    if vcf_samples != panel_samples:
        only_vcf = sorted(vcf_samples - panel_samples)
        only_panel = sorted(panel_samples - vcf_samples)
        raise RuntimeError(f"VCF and panel samples differ: {only_vcf=}, {only_panel=}")

    # Split target and reference samples.
    print("Splitting samples into target and ref.")
    random.seed(args.seed)
    target_samples = set()
    for pop, samples in pops.items():
        target_samples.update(random.sample(samples, args.targets_per_pop))
    ref_samples = panel_samples - target_samples

    args.output_dir.mkdir(parents=True, exist_ok=True)

    ref_list_file = args.output_dir / "ref_samples.txt"
    ref_list_file.write_text("\n".join(sorted(ref_samples)) + "\n")

    target_list_file = args.output_dir / "target_samples.txt"
    target_list_file.write_text("\n".join(sorted(target_samples)) + "\n")

    # Split vcf.
    ref_vcf = args.output_dir / "ref.vcf.gz"
    target_vcf = args.output_dir / "target.vcf.gz"
    with tempfile.TemporaryDirectory(prefix="split-", dir=args.output_dir) as tmp:
        ref_all = f"{tmp}/ref.all.bcf"
        target_all = f"{tmp}/target.all.bcf"

        # Subset first so that AC/AN and the minor-allele filter describe the
        # reference panel, rather than all 2,504 individuals.
        # -S keeps only samples in the newline separated sample name file.
        print("Getting ref sequences.")
        subprocess.run([
            "bcftools", "view",
            "--threads", str(args.threads),
            "-S", str(ref_list_file),
            "-Ob", "-o", str(ref_all),
            str(args.vcf),
        ], check=True)

        # -c 1:minor sets the minimum minor allele count to 1.
        print("Removing monomorphic ref records.")
        subprocess.run([
            "bcftools", "view",
            "--threads", str(args.threads),
            "-c", "1:minor",
            "-Oz", "-o", str(ref_vcf),
            str(ref_all),
        ], check=True)
        subprocess.run(["bcftools", "index", "--force", str(ref_vcf)], check=True)

        print("Getting target sequences.")
        subprocess.run([
            "bcftools", "view",
            "--threads", str(args.threads),
            "-S", str(target_list_file),
            "-Ob", "-o", str(target_all),
            str(args.vcf),
        ], check=True)
        subprocess.run(["bcftools", "index", "--force", str(target_all)], check=True)

        # isec intersects variant loc rows, keeps only variants that -n=2 files have,
        # and (-w1) write records from 1, target all (1 based idx).
        print("Matching target records to polymorphic ref records.")
        subprocess.run([
            "bcftools", "isec",
            "--threads", str(args.threads),
            "-n=2",
            "-w1",
            "-Oz", "-o", str(target_vcf),
            str(target_all),
            str(ref_vcf)
        ], check=True)
        subprocess.run(["bcftools", "index", "--force", str(target_vcf)], check=True)

    print("Checking # target and ref records.")
    ref_records = record_count(ref_vcf)
    target_records = record_count(target_vcf)
    if ref_records != target_records:
        raise RuntimeError(f"ref and target record counts differ: {ref_records} != {target_records}")
