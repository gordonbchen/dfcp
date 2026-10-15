#!/usr/bin/env python3

import argparse
import os
import signal
import subprocess
from dataclasses import dataclass
from pathlib import Path


DATA = Path("data/1000g_phase3_v5b")
REF = DATA / "ref_target/ref.vcf.gz"
OBS = DATA / "mask_target/target_observed.vcf.gz"
TRUE = DATA / "mask_target/target_masked_true.vcf.gz"
MAP = DATA / "plink.chr20.GRCh37.map"
OUT = DATA / "windows/beagle"
JAR = Path("beagle/beagle.27Feb25.75f.jar")


@dataclass
class Stats:
    n_loci: int = 0
    n: int = 0
    n_correct: int = 0
    sum_q: float = 0.0
    sum_y: int = 0
    sum_qq: float = 0.0
    sum_qy: float = 0.0

    def add(self, q: float, y: int) -> None:
        self.n += 1
        self.n_correct += (q >= 0.5) == y
        self.sum_q += q
        self.sum_y += y
        self.sum_qq += q * q
        self.sum_qy += q * y

    def r2(self) -> float:
        cov = self.n * self.sum_qy - self.sum_q * self.sum_y
        var_q = self.n * self.sum_qq - self.sum_q * self.sum_q
        var_y = self.n * self.sum_y - self.sum_y * self.sum_y
        if var_q <= 0.0 or var_y <= 0:
            return -1.0
        return min(1.0, max(0.0, cov * cov / (var_q * var_y)))


def query(vcf: Path, fmt: str, include: str | None = None) -> subprocess.Popen:
    command = ["bcftools", "query", "-f", fmt]
    if include:
        command += ["-i", include]
    return subprocess.Popen(command + [vcf], stdout=subprocess.PIPE, text=True)


def rows(process: subprocess.Popen):
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip().split("\t")
        yield (fields[0], int(fields[1]), fields[2], fields[3]), fields


def advance(row, iterator, key):
    while row[0] < key:
        row = next(iterator)
    return row


def evaluate(vcf: Path, output: Path) -> None:
    ref_process = query(REF, "%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n")
    true_process = query(TRUE, "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n")
    prob_process = query(vcf, "%CHROM\t%POS\t%REF\t%ALT[\t%AP1\t%AP2]\n", "INFO/IMP=1")
    processes = [ref_process, true_process, prob_process]
    ref_rows, true_rows, prob_rows = map(rows, processes)
    stats: dict[int, Stats] = {}

    try:
        ref, truth, prob = next(ref_rows), next(true_rows), next(prob_rows)
        while True:
            key = max(ref[0], truth[0], prob[0])
            ref = advance(ref, ref_rows, key)
            truth = advance(truth, true_rows, key)
            prob = advance(prob, prob_rows, key)
            if ref[0] != truth[0] or ref[0] != prob[0]:
                continue

            ac, an = map(int, ref[1][4:])
            stat = stats.setdefault(min(ac, an - ac), Stats())
            stat.n_loci += 1
            alt_is_minor = ac <= an - ac
            probabilities = iter(map(float, prob[1][4:]))
            for gt in truth[1][4:]:
                q1, q2 = next(probabilities), next(probabilities)
                for q, y in zip((q1, q2), map(int, gt.split("|"))):
                    stat.add(q if alt_is_minor else 1.0 - q, y if alt_is_minor else 1 - y)
            ref, truth, prob = next(ref_rows), next(true_rows), next(prob_rows)
    except StopIteration:
        pass
    finally:
        for process in processes:
            if process.stdout:
                process.stdout.close()
            process.wait()
            if process.returncode not in (0, -signal.SIGPIPE):
                raise subprocess.CalledProcessError(process.returncode, process.args)

    with output.open("w") as stream:
        stream.write("mac\tn_loci\tn_predictions\tr2\taccuracy\n")
        for mac, stat in sorted(stats.items()):
            stream.write(
                f"{mac}\t{stat.n_loci}\t{stat.n}\t{stat.r2():.10g}"
                f"\t{stat.n_correct / stat.n:.10g}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run and evaluate full-chromosome Beagle imputation.")
    parser.add_argument("--nthreads", type=int, default=os.cpu_count())
    parser.add_argument("--heap", default="16g")
    parser.add_argument("--eval-only", action="store_true")
    args = parser.parse_args()

    if not args.eval_only:
        subprocess.run([
            "java", f"-Xmx{args.heap}", "-jar", str(JAR), f"gt={OBS}", f"ref={REF}",
            f"map={MAP}", f"out={OUT}", "ap=true", "seed=1", f"nthreads={args.nthreads}",
        ], check=True)
    output = OUT.parent / "eval_impute_beagle.tsv"
    evaluate(OUT.with_suffix(".vcf.gz"), output)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
