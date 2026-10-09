#!/usr/bin/env python3

"""Plot observed-allele agreement between each target and reference haplotype."""

import argparse
import struct
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from plotly_html import ensure_plotly_asset


def read_sequences(path: Path) -> tuple[int, int, np.ndarray]:
    data = path.read_bytes()
    if data[:4] != b"DFCP":
        raise ValueError(f"invalid sequence file: {path}")
    n_sequences, n_loci = struct.unpack_from("<II", data, 4)
    words_per_locus = (n_sequences + 63) // 64
    words = np.frombuffer(data, dtype="<u8", offset=12)
    if words.size != n_loci * words_per_locus:
        raise ValueError(f"wrong sequence file size: {path}")
    return n_sequences, n_loci, words.reshape(n_loci, words_per_locus)


def pack_sequences(words: np.ndarray, n_sequences: int) -> np.ndarray:
    # The file is locus-major; each output row is one haplotype over observed loci.
    alleles = np.unpackbits(words.view(np.uint8), axis=1, bitorder="little")[:, :n_sequences]
    return np.packbits(alleles.T, axis=1, bitorder="little")


def window_match_counts(window: Path) -> tuple[int, np.ndarray]:
    n_ref, n_loci, ref_words = read_sequences(window / "ref.bin")
    n_target, n_observed, target_words = read_sequences(window / "target_observed.bin")
    observed = np.loadtxt(window / "observed_loci.txt", dtype=np.int64, ndmin=1)
    if observed.size != n_observed or np.any(observed < 0) or np.any(observed >= n_loci):
        raise ValueError(f"invalid observed loci: {window}")

    ref = pack_sequences(ref_words[observed], n_ref)
    target = pack_sequences(target_words, n_target)
    counts = np.empty((n_target, n_ref), dtype=np.uint16)
    for start in range(0, n_target, 8):
        end = min(start + 8, n_target)
        mismatches = np.bitwise_count(target[start:end, None, :] ^ ref[None, :, :])
        counts[start:end] = n_observed - mismatches.sum(axis=2, dtype=np.uint16)
    return n_observed, counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("windows", type=Path)
    parser.add_argument("--output", type=Path, default=Path("output/ref_match_hist.html"))
    args = parser.parse_args()

    windows = sorted(args.windows.glob("window_*"))
    if not windows:
        raise ValueError(f"no windows found in {args.windows}")

    count_hist = np.zeros(1, dtype=np.int64)
    percent_hist = np.zeros(101, dtype=np.int64)
    retained = {fraction: [] for fraction in (0.70, 0.75, 0.80, 0.85, 0.90, 0.95)}
    n_observed_range = [float("inf"), 0]
    for window in windows:
        n_observed, counts = window_match_counts(window)
        n_observed_range[0] = min(n_observed_range[0], n_observed)
        n_observed_range[1] = max(n_observed_range[1], n_observed)
        histogram = np.bincount(counts.ravel(), minlength=n_observed + 1)
        if histogram.size > count_hist.size:
            count_hist.resize(histogram.size, refcheck=False)
        count_hist[:histogram.size] += histogram
        percentages = 100 * counts.astype(np.float32) / n_observed
        percent_hist += np.bincount(np.minimum(percentages.astype(np.int16), 100).ravel(),
                                    minlength=101)
        for fraction in retained:
            retained[fraction].extend(np.count_nonzero(counts >= fraction * n_observed, axis=1))
        print(f"{window.name}: {n_observed} observed sites, {counts.shape[0]} targets")

    n_pairs = int(count_hist.sum())
    n_targets = len(next(iter(retained.values())))
    figure = make_subplots(rows=2, cols=1, vertical_spacing=0.14,
                           subplot_titles=("Matching observed alleles (raw count)",
                                           "Matching observed alleles (percent)"))
    figure.add_trace(go.Bar(x=np.flatnonzero(count_hist), y=count_hist[count_hist > 0],
                            marker_color="#4063a0", name="Raw count"), row=1, col=1)
    figure.add_trace(go.Bar(x=np.arange(101), y=percent_hist, marker_color="#b46a3c",
                            name="Percent"), row=2, col=1)
    figure.update_layout(
        title=f"Target–reference agreement: {len(windows)} windows, {n_targets:,} targets, "
              f"{n_pairs:,} pairs; {n_observed_range[0]}–{n_observed_range[1]} observed sites/window",
        template="plotly_white", showlegend=False, height=850, bargap=0,
        margin={"l": 90, "r": 40, "t": 115, "b": 70},
    )
    figure.update_xaxes(title_text="Matching sites", row=1, col=1)
    figure.update_xaxes(title_text="Matching sites (%)", row=2, col=1)
    figure.update_yaxes(title_text="Target–reference pairs", row=1, col=1)
    figure.update_yaxes(title_text="Target–reference pairs", row=2, col=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(args.output, include_plotlyjs=ensure_plotly_asset(args.output))

    print(f"Wrote {args.output}")
    for fraction, values in retained.items():
        values = np.asarray(values)
        print(f">= {fraction:.0%} matching: median {np.median(values):.0f}, "
              f"p10 {np.percentile(values, 10):.0f}, p90 {np.percentile(values, 90):.0f}, "
              f"zero-targets {np.count_nonzero(values == 0)}")


if __name__ == "__main__":
    main()
