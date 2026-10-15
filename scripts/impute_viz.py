#!/usr/bin/env python3

import argparse
import csv
import math
import sys
from pathlib import Path

import plotly.graph_objects as go
from plotly.colors import qualitative
from plotly_html import ensure_plotly_asset

FIELDS = ["mac", "n_loci", "n_predictions", "r2", "accuracy"]
Evaluation = tuple[list[int], list[int], list[float], list[float]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot imputation accuracy and r-squared by reference minor-allele count.",
    )
    parser.add_argument("evaluations", type=Path, nargs="+", help="aggregate imputation TSV files")
    parser.add_argument("--output", type=Path, default=Path("impute.html"))
    return parser.parse_args()


def read_evaluation(path: Path) -> Evaluation:
    macs = []
    n_loci = []
    r2 = []
    accuracy = []
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if rows.fieldnames != FIELDS:
            raise ValueError(f"unexpected imputation evaluation header: {path}")
        for row_number, row in enumerate(rows, start=2):
            try:
                values = (
                    int(row["mac"]),
                    int(row["n_loci"]),
                    float(row["r2"]),
                    float(row["accuracy"]),
                )
            except ValueError as error:
                raise ValueError(f"invalid value on {path} row {row_number}") from error
            mac, loci, row_r2, row_accuracy = values
            if mac < 0 or (macs and mac <= macs[-1]) or loci <= 0:
                raise ValueError(f"invalid counts on {path} row {row_number}")
            if not math.isfinite(row_r2) or (row_r2 < 0.0 and row_r2 != -1.0) or row_r2 > 1.0:
                raise ValueError(f"invalid r-squared on {path} row {row_number}")
            if not math.isfinite(row_accuracy) or not 0.0 <= row_accuracy <= 1.0:
                raise ValueError(f"invalid accuracy on {path} row {row_number}")
            macs.append(mac)
            n_loci.append(loci)
            r2.append(row_r2)
            accuracy.append(row_accuracy)
    if not macs:
        raise ValueError(f"imputation evaluation contains no MAC bins: {path}")
    return macs, n_loci, r2, accuracy


def make_figure(evaluations: list[tuple[Path, Evaluation]]) -> go.Figure:
    hover = (
        "reference MAC=%{x}<br>%{fullData.name}=%{y:.4f}"
        "<br>n_loci=%{customdata:,}<extra></extra>"
    )
    figure = go.Figure()
    for i, (path, (macs, n_loci, r2, accuracy)) in enumerate(evaluations):
        color = qualitative.Plotly[i % len(qualitative.Plotly)]
        for metric, values, dash in (("r²", r2, "solid"), ("accuracy", accuracy, "dash")):
            figure.add_trace(go.Scatter(
                x=macs,
                y=[None if value < 0.0 else value for value in values],
                customdata=n_loci,
                mode="lines",
                name=f"{path.name} {metric}",
                line={"color": color, "width": 2, "dash": dash},
                hovertemplate=hover,
            ))
    figure.update_layout(
        title="Imputation performance by reference minor-allele count",
        xaxis_title="Reference minor-allele count (haplotypes)",
        xaxis_type="log",
        yaxis_title="Imputation metric",
        yaxis_range=[0, 1],
        template="plotly_white",
        hovermode="x unified",
        legend={"y": 1, "x": 1.02, "xanchor": "left"},
        margin={"l": 70, "r": 300, "t": 80, "b": 65},
    )
    return figure


def main() -> None:
    args = parse_args()
    evaluations = [(path, read_evaluation(path)) for path in args.evaluations]
    figure = make_figure(evaluations)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(
        args.output,
        include_plotlyjs=ensure_plotly_asset(args.output),
        config={"responsive": True, "displaylogo": False},
    )
    print(
        f"wrote {args.output}: evaluations={len(evaluations)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
