#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def write_vcf(gen_file: Path, output: Path) -> None:
    with gen_file.open() as source, output.open("w") as vcf:
        header = source.readline().split()
        if header[:4] != ["Chrom", "Pos", "Anc_all", "Der_all"] or (len(header) - 4) % 2:
            raise ValueError("expected an even number of haplotypes in a fastsimcoal .gen file")

        samples = [f"sample_{i}" for i in range((len(header) - 4) // 2)]
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        previous_pos = 0
        for i, line in enumerate(source):
            fields = line.split()
            if not fields:
                continue
            chrom, pos, ref, alt = fields[:4]
            haps = fields[4:]
            if len(fields) != len(header) or any(allele not in ("0", "1") for allele in haps):
                raise ValueError(f"invalid genotype row {i + 2}")
            beagle_pos = max(int(pos), previous_pos + 1)
            previous_pos = beagle_pos
            genotypes = ["|".join(haps[j:j + 2]) for j in range(0, len(haps), 2)]
            vcf.write(f"{chrom}\t{beagle_pos}\tfsc_{i}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t")
            vcf.write("\t".join(genotypes) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract Beagle 4 DAG-edge R assignments for a fastsimcoal run."
    )
    parser.add_argument("gen_file", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--source", type=Path, default=ROOT / "beagle/src")
    parser.add_argument("--javac", default="javac")
    parser.add_argument("--java", default="java")
    args = parser.parse_args()

    rassign = args.output_dir / "beagle4.r_assign.bin"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sources = sorted(args.source.rglob("*.java"))

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        classes = temp_dir / "classes"
        vcf = temp_dir / "input.vcf"
        subprocess.run([args.javac, "-d", str(classes), *map(str, sources)], check=True)
        write_vcf(args.gen_file, vcf)
        subprocess.run([
            args.java, "-cp", str(classes), "main.Main", f"gt={vcf}",
            f"out={args.output_dir / 'beagle4'}", "impute=false", "gprobs=false", "usephase=true",
            "burnin-its=0", "phase-its=0", "impute-its=0", "window=50000", "overlap=0",
            f"rassign={rassign}",
        ], check=True)

    print(rassign)


if __name__ == "__main__":
    main()
