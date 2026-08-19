#!/usr/bin/env python3

import argparse
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
import struct
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
PREP_DIR = REPO_ROOT / "data/1000g_phase3_v5b/dfcp_prep"
DEFAULT_REF = PREP_DIR / "ref.bin"
DEFAULT_TARGET_OBSERVED = PREP_DIR / "target_observed.bin"
DEFAULT_TARGET_MASKED_TRUE = PREP_DIR / "target_masked_true.bin"
DEFAULT_OBSERVED_LOCI = PREP_DIR / "observed_loci.txt"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data/1000g_phase3_v5b/windows"

SEQ_HEADER = struct.Struct("<4sII")
SEQ_MAGIC = b"DFCP"
COPY_BUFFER_SIZE = 8 * 1024 * 1024


@dataclass(frozen=True)
class SeqFile:
    path: Path
    n_sequences: int
    n_loci: int
    bytes_per_locus: int


@dataclass(frozen=True)
class Window:
    index: int
    start: int
    end: int

    @property
    def size(self) -> int:
        return self.end - self.start


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


def read_seq_file(path: Path) -> SeqFile:
    with path.open("rb") as stream:
        header = stream.read(SEQ_HEADER.size)
    if len(header) != SEQ_HEADER.size:
        raise ValueError(f"truncated sequence header: {path}")

    magic, n_sequences, n_loci = SEQ_HEADER.unpack(header)
    if magic != SEQ_MAGIC:
        raise ValueError(f"invalid sequence magic: {path}")
    if n_sequences == 0 or n_loci == 0:
        raise ValueError(f"sequence dimensions must be positive: {path}")

    bytes_per_locus = ((n_sequences + 63) // 64) * 8
    expected_size = SEQ_HEADER.size + n_loci * bytes_per_locus
    actual_size = path.stat().st_size
    if actual_size != expected_size:
        raise ValueError(f"sequence file has {actual_size} bytes; expected {expected_size}: {path}")
    return SeqFile(path, n_sequences, n_loci, bytes_per_locus)


def read_observed_loci(path: Path, n_loci: int) -> list[int]:
    loci = []
    with path.open() as stream:
        for row, line in enumerate(stream, start=1):
            value = line.strip()
            if not value:
                raise ValueError(f"blank observed-locus row {row}: {path}")
            try:
                locus = int(value)
            except ValueError as error:
                raise ValueError(f"invalid observed locus on row {row}: {value}") from error
            if not 0 <= locus < n_loci:
                raise ValueError(f"observed locus {locus} on row {row} is outside [0, {n_loci})")
            if loci and locus <= loci[-1]:
                raise ValueError("observed loci must be unique and strictly increasing")
            loci.append(locus)
    return loci


def copy_exact(source, destination, n_bytes: int) -> None:
    remaining = n_bytes
    while remaining:
        block = source.read(min(remaining, COPY_BUFFER_SIZE))
        if not block:
            raise ValueError("sequence file ended while copying a validated locus range")
        destination.write(block)
        remaining -= len(block)


def write_seq_slice(seq_file: SeqFile, start: int, end: int, output: Path) -> None:
    if not 0 <= start < end <= seq_file.n_loci:
        raise ValueError(f"invalid locus range [{start}, {end}) for {seq_file.path}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with seq_file.path.open("rb") as source, output.open("wb") as destination:
        destination.write(SEQ_HEADER.pack(SEQ_MAGIC, seq_file.n_sequences, end - start))
        source.seek(SEQ_HEADER.size + start * seq_file.bytes_per_locus)
        copy_exact(source, destination, (end - start) * seq_file.bytes_per_locus)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split bitpacked DFCP chromosome inputs into equal-locus overlapping windows.",
    )
    parser.add_argument("--ref", type=Path, default=DEFAULT_REF)
    parser.add_argument("--target-observed", type=Path, default=DEFAULT_TARGET_OBSERVED)
    parser.add_argument("--target-masked-true", type=Path, default=DEFAULT_TARGET_MASKED_TRUE)
    parser.add_argument("--observed-loci", type=Path, default=DEFAULT_OBSERVED_LOCI)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-windows", type=int, required=True)
    parser.add_argument("--overlap", type=int, default=0, help="reference loci shared by adjacent windows")
    parser.add_argument("--first-only", action="store_true", help="write only window 0")
    return parser.parse_args()


def validate_inputs(
    ref: SeqFile,
    target_observed: SeqFile,
    target_masked_true: SeqFile,
    observed_loci: list[int],
) -> None:
    if target_observed.n_sequences != target_masked_true.n_sequences:
        raise ValueError("observed target and masked truth have different sequence counts")
    if target_observed.n_loci != len(observed_loci):
        raise ValueError("observed-loci rows do not match the observed-target locus count")
    if target_masked_true.n_loci != ref.n_loci - len(observed_loci):
        raise ValueError("masked-truth loci are not the complement of the observed loci")


def require_empty_output_dir(path: Path) -> None:
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"output directory must not exist or must be empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def write_window(
    output_dir: Path,
    window: Window,
    ref: SeqFile,
    target_observed: SeqFile,
    target_masked_true: SeqFile,
    observed_loci: list[int],
    overlap: int,
) -> tuple[int, int, int]:
    first_observed = bisect_left(observed_loci, window.start)
    after_observed = bisect_left(observed_loci, window.end)
    observed_count = after_observed - first_observed
    first_masked = window.start - first_observed
    after_masked = window.end - after_observed
    masked_count = after_masked - first_masked
    if not observed_count:
        raise ValueError(f"window {window.index} has no observed target loci")
    if not masked_count:
        raise ValueError(f"window {window.index} has no masked truth loci")

    window_dir = output_dir / f"window_{window.index:04d}"
    window_dir.mkdir()
    write_seq_slice(ref, window.start, window.end, window_dir / "ref.bin")
    write_seq_slice(
        target_observed, first_observed, after_observed, window_dir / "target_observed.bin"
    )
    write_seq_slice(
        target_masked_true, first_masked, after_masked, window_dir / "target_masked_true.bin"
    )

    local_observed = observed_loci[first_observed:after_observed]
    observed_path = window_dir / "observed_loci.txt"
    observed_path.write_text("".join(f"{locus - window.start}\n" for locus in local_observed))

    output_bytes = sum(path.stat().st_size for path in window_dir.iterdir())
    overlap_prev = overlap if window.index else 0
    print(
        f"window={window.index:04d} ref=[{window.start},{window.end}) loci={window.size} "
        f"observed={observed_count} masked={masked_count} overlap_prev={overlap_prev} "
        f"bytes={output_bytes}",
        file=sys.stderr,
    )
    return observed_count, masked_count, output_bytes


def main() -> None:
    args = parse_args()
    ref = read_seq_file(args.ref)
    target_observed = read_seq_file(args.target_observed)
    target_masked_true = read_seq_file(args.target_masked_true)
    observed_loci = read_observed_loci(args.observed_loci, ref.n_loci)
    validate_inputs(ref, target_observed, target_masked_true, observed_loci)
    windows = make_windows(ref.n_loci, args.n_windows, args.overlap)
    selected = windows[:1] if args.first_only else windows
    require_empty_output_dir(args.output_dir)

    manifest_path = args.output_dir / "windows.tsv"
    total_bytes = 0
    with manifest_path.open("w") as manifest:
        manifest.write("window\tstart\tend\tloci\tobserved\tmasked\toverlap_previous\n")
        for window in selected:
            observed, masked, output_bytes = write_window(
                args.output_dir,
                window,
                ref,
                target_observed,
                target_masked_true,
                observed_loci,
                args.overlap,
            )
            overlap_prev = args.overlap if window.index else 0
            manifest.write(
                f"{window.index}\t{window.start}\t{window.end}\t{window.size}\t{observed}\t"
                f"{masked}\t{overlap_prev}\n"
            )
            total_bytes += output_bytes

    print(
        f"wrote {len(selected)} of {len(windows)} windows to {args.output_dir}; "
        f"window_bytes={total_bytes} manifest={manifest_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
