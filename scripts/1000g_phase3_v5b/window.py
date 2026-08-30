#!/usr/bin/env python3

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DATA_DIR = REPO_ROOT / "data/1000g_phase3_v5b"
DEFAULT_REF = DATA_DIR / "ref_target/ref.vcf.gz"
DEFAULT_TARGET_OBSERVED = DATA_DIR / "mask_target/target_observed.vcf.gz"
DEFAULT_TARGET_MASKED_TRUE = DATA_DIR / "mask_target/target_masked_true.vcf.gz"
DEFAULT_OUTPUT_DIR = DATA_DIR / "windows"


@dataclass(frozen=True)
class Window:
    index: int
    start: int
    end: int

    @property
    def size(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Position:
    chromosome: str
    value: int


def make_windows(n_loci: int, n_windows: int, overlap: int) -> list[Window]:
    if n_loci <= 0:
        raise ValueError("the number of loci must be positive")
    if n_windows <= 0:
        raise ValueError("the number of windows must be positive")
    if overlap < 0:
        raise ValueError("overlap must be nonnegative")

    total_memberships = n_loci + (n_windows - 1) * overlap
    base_size, remainder = divmod(total_memberships, n_windows)
    if n_windows > 1 and base_size <= overlap:
        raise ValueError("the overlap leaves no new loci in at least one window")

    windows = []
    start = 0
    for index in range(n_windows):
        size = base_size + (index < remainder)
        end = start + size
        windows.append(Window(index, start, end))
        start = end - overlap

    if windows[-1].end != n_loci:
        raise RuntimeError("internal window construction error")
    return windows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split aligned reference and target VCFs into equal-locus overlapping windows.",
    )
    parser.add_argument("--ref", type=Path, default=DEFAULT_REF)
    parser.add_argument("--target-observed", type=Path, default=DEFAULT_TARGET_OBSERVED)
    parser.add_argument("--target-masked-true", type=Path, default=DEFAULT_TARGET_MASKED_TRUE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-windows", type=int, required=True)
    parser.add_argument("--overlap", type=int, default=0, help="reference loci shared by neighbors")
    parser.add_argument(
        "--n-generate",
        type=int,
        help="generate only the first N windows; the default generates every window",
    )
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.n_generate is None:
        args.n_generate = args.n_windows
    if not 1 <= args.n_generate <= args.n_windows:
        parser.error("--n-generate must be between 1 and --n-windows")
    if args.threads <= 0:
        parser.error("--threads must be positive")
    return args


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def record_count(path: Path) -> int:
    result = subprocess.run(
        ["bcftools", "index", "--nrecords", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout)


def read_boundary_positions(vcf: Path, boundary_indexes: list[int]) -> dict[int, Position]:
    process = subprocess.Popen(
        ["bcftools", "query", "-f", "%CHROM\t%POS\n", str(vcf)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    wanted = set(boundary_indexes)
    positions = {}
    for index, line in enumerate(process.stdout):
        if index in wanted:
            chromosome, position = line.rstrip().split("\t")
            positions[index] = Position(chromosome, int(position))
    if process.wait():
        raise subprocess.CalledProcessError(process.returncode, process.args)
    if len(positions) != len(wanted):
        raise RuntimeError("failed to find every window boundary in the reference VCF")
    return positions


def write_region(source: Path, output: Path, region: str, threads: int) -> int:
    run([
        "bcftools", "view",
        "--threads", str(threads),
        "--regions", region,
        "-Oz", "-o", str(output),
        str(source),
    ])
    run(["bcftools", "index", "--force", "--threads", str(threads), str(output)])
    return record_count(output)


def require_empty_output_dir(path: Path) -> None:
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"output directory must not exist or must be empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def write_window(
    args: argparse.Namespace,
    window: Window,
    first: Position,
    last: Position,
) -> tuple[int, int]:
    if first.chromosome != last.chromosome:
        raise ValueError("a window cannot cross chromosomes")
    region = f"{first.chromosome}:{first.value}-{last.value}"
    window_dir = args.output_dir / f"window_{window.index:04d}"
    window_dir.mkdir()
    n_ref = write_region(args.ref, window_dir / "ref.vcf.gz", region, args.threads)
    observed = write_region(
        args.target_observed, window_dir / "target_observed.vcf.gz", region, args.threads
    )
    masked = write_region(
        args.target_masked_true, window_dir / "target_masked_true.vcf.gz", region, args.threads
    )
    if n_ref != window.size or observed == 0 or masked == 0 or observed + masked != n_ref:
        raise RuntimeError(
            f"window {window.index} expected {window.size} records but has "
            f"ref={n_ref}, observed={observed}, masked={masked}"
        )
    print(
        f"window={window.index:04d} ref=[{window.start},{window.end}) loci={n_ref} "
        f"observed={observed} masked={masked}",
        file=sys.stderr,
    )
    return observed, masked


def main() -> None:
    args = parse_args()
    n_loci = record_count(args.ref)
    windows = make_windows(n_loci, args.n_windows, args.overlap)
    boundary_indexes = [index for window in windows for index in (window.start, window.end - 1)]
    boundaries = read_boundary_positions(args.ref, boundary_indexes)
    require_empty_output_dir(args.output_dir)

    generated: dict[int, tuple[int, int]] = {}
    for window in windows[:args.n_generate]:
        generated[window.index] = write_window(
            args, window, boundaries[window.start], boundaries[window.end - 1]
        )

    manifest = args.output_dir / "windows.tsv"
    with manifest.open("w") as stream:
        stream.write(
            "window\tstart\tend\tplanned_loci\tchrom\tfirst_pos\tlast_pos\tgenerated"
            "\tobserved\tmasked\toverlap_previous\n"
        )
        for window in windows:
            first = boundaries[window.start]
            last = boundaries[window.end - 1]
            counts = generated.get(window.index)
            observed, masked = counts if counts else (".", ".")
            overlap = args.overlap if window.index else 0
            stream.write(
                f"{window.index}\t{window.start}\t{window.end}\t{window.size}\t"
                f"{first.chromosome}\t{first.value}\t{last.value}\t{int(counts is not None)}\t"
                f"{observed}\t{masked}\t{overlap}\n"
            )

    print(
        f"wrote {args.n_generate} of {args.n_windows} VCF windows to {args.output_dir}; "
        f"manifest={manifest}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
