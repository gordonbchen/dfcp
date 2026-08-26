#!/usr/bin/env python3

import argparse
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly_html import ensure_plotly_asset

EVAL_HEADER = struct.Struct("<4sI")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot per-locus imputation r-squared against reference minor-allele count.",
    )
    parser.add_argument("evaluation", type=Path)
    parser.add_argument("reference", type=Path, help="windowed reference VCF")
    parser.add_argument("observed_loci", type=Path)
    parser.add_argument("--output", type=Path, default=Path("impute.html"))
    parser.add_argument("--max-points", type=int, default=200_000)
    args = parser.parse_args()
    if args.max_points <= 0:
        parser.error("--max-points must be positive")
    return args


def read_reference_mac(path: Path) -> np.ndarray:
    command = ["bcftools", "query", "-f", "%AC\\t%AN\\n", str(path)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    counts = []
    for row, line in enumerate(result.stdout.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) != 2 or "," in fields[0]:
            raise ValueError(f"reference VCF row {row} is not biallelic: {line}")
        ac, an = map(int, fields)
        if not 0 <= ac <= an:
            raise ValueError(f"invalid AC/AN on reference VCF row {row}: {ac}/{an}")
        counts.append(min(ac, an - ac))
    if not counts:
        raise ValueError(f"reference VCF contains no variants: {path}")
    return np.asarray(counts, dtype=np.int32)


def read_observed_loci(path: Path, n_loci: int) -> np.ndarray:
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
                raise ValueError(f"observed locus {locus} is outside [0, {n_loci})")
            if loci and locus <= loci[-1]:
                raise ValueError("observed loci must be unique and strictly increasing")
            loci.append(locus)
    return np.asarray(loci, dtype=np.int64)


def read_evaluation(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open("rb") as stream:
        header = stream.read(EVAL_HEADER.size)
    if len(header) != EVAL_HEADER.size:
        raise ValueError(f"truncated imputation evaluation header: {path}")
    magic, n_loci = EVAL_HEADER.unpack(header)
    if magic != b"DFIE" or n_loci == 0:
        raise ValueError(f"invalid imputation evaluation header: {path}")
    expected_size = EVAL_HEADER.size + 2 * n_loci * 4
    if path.stat().st_size != expected_size:
        raise ValueError(f"invalid imputation evaluation file size: {path}")
    values = np.memmap(path, dtype="<f4", mode="r", offset=EVAL_HEADER.size, shape=(2, n_loci))
    if not np.all(np.isfinite(values)):
        raise ValueError(f"imputation evaluation contains non-finite values: {path}")
    if np.any((values[0] < 0.0) & (values[0] != -1.0)) or np.any(values[0] > 1.0):
        raise ValueError(f"imputation evaluation contains invalid r-squared values: {path}")
    if np.any((values[1] < 0.0) | (values[1] > 1.0)):
        raise ValueError(f"imputation evaluation contains invalid accuracy values: {path}")
    return values[0], values[1]


def masked_indexes(n_loci: int, observed: np.ndarray, n_masked: int) -> np.ndarray:
    is_observed = np.zeros(n_loci, dtype=bool)
    is_observed[observed] = True
    masked = np.flatnonzero(~is_observed)
    if masked.size != n_masked:
        raise ValueError(
            f"evaluation has {n_masked} loci but the observed-loci complement has {masked.size}"
        )
    return masked


def make_figure(
    r2: np.ndarray,
    accuracy: np.ndarray,
    minor_counts: np.ndarray,
    max_points: int,
) -> tuple[go.Figure, int]:
    valid = np.flatnonzero(r2 >= 0.0)
    if valid.size > max_points:
        rng = np.random.default_rng(0)
        plotted = np.sort(rng.choice(valid, max_points, replace=False))
    else:
        plotted = valid

    mac_totals = np.bincount(minor_counts[valid], weights=r2[valid])
    mac_loci = np.bincount(minor_counts[valid])
    populated = np.flatnonzero(mac_loci)
    mac_means = mac_totals[populated] / mac_loci[populated]

    figure = go.Figure()
    figure.add_trace(go.Scattergl(
        x=minor_counts[plotted],
        y=r2[plotted],
        mode="markers",
        name="Masked loci",
        marker={"color": "#3366cc", "opacity": 0.22, "size": 4},
        hovertemplate="reference MAC=%{x}<br>r²=%{y:.4f}<extra></extra>",
    ))
    figure.add_trace(go.Scatter(
        x=populated,
        y=mac_means,
        mode="lines",
        name="Mean at each MAC",
        line={"color": "#ef8a35", "width": 2.5},
        hovertemplate="reference MAC=%{x}<br>mean r²=%{y:.4f}<extra></extra>",
    ))
    figure.update_layout(
        title="Imputation r² by reference minor-allele count",
        xaxis_title="Reference minor-allele count (haplotypes)",
        yaxis_title="Imputation r² across target haplotypes",
        template="plotly_white",
        hovermode="closest",
        legend={"orientation": "h", "y": 1.08, "x": 1, "xanchor": "right"},
        margin={"l": 70, "r": 30, "t": 80, "b": 65},
    )
    mean_r2 = float(np.mean(r2[valid])) if valid.size else -1.0
    figure.add_annotation(
        x=0,
        y=1.08,
        xref="paper",
        yref="paper",
        showarrow=False,
        text=(
            f"mean r²={mean_r2:.4f}; mean accuracy={float(np.mean(accuracy)):.4f}; "
            f"undefined r²={r2.size - valid.size:,}"
        ),
    )
    return figure, plotted.size


def main() -> None:
    args = parse_args()
    reference_mac = read_reference_mac(args.reference)
    n_loci = reference_mac.size
    r2, accuracy = read_evaluation(args.evaluation)
    observed = read_observed_loci(args.observed_loci, n_loci)
    masked = masked_indexes(n_loci, observed, r2.size)
    minor_counts = reference_mac[masked]
    figure, n_plotted = make_figure(r2, accuracy, minor_counts, args.max_points)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(
        args.output,
        include_plotlyjs=ensure_plotly_asset(args.output),
        config={"responsive": True, "displaylogo": False},
    )
    n_defined = int(np.count_nonzero(r2 >= 0.0))
    print(
        f"wrote {args.output}: loci={r2.size} defined_r2={n_defined} plotted={n_plotted} "
        f"mean_r2={float(np.mean(r2[r2 >= 0.0])) if n_defined else -1.0:.6f}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
