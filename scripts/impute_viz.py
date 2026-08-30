#!/usr/bin/env python3

import argparse
import csv
import math
import sys
from pathlib import Path

import plotly.graph_objects as go
from plotly_html import ensure_plotly_asset

FIELDS = ["mac", "n_loci", "n_predictions", "r2", "accuracy"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot pooled imputation accuracy and r-squared by reference minor-allele count.",
    )
    parser.add_argument("evaluation", type=Path, help="aggregate TSV written by eval_impute")
    parser.add_argument("--output", type=Path, default=Path("impute.html"))
    return parser.parse_args()


def read_evaluation(path: Path) -> tuple[list[int], list[int], list[int], list[float], list[float]]:
    macs = []
    n_loci = []
    n_predictions = []
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
                    int(row["n_predictions"]),
                    float(row["r2"]),
                    float(row["accuracy"]),
                )
            except ValueError as error:
                raise ValueError(f"invalid value on {path} row {row_number}") from error
            mac, loci, predictions, row_r2, row_accuracy = values
            if mac < 0 or (macs and mac <= macs[-1]) or loci <= 0 or predictions <= 0:
                raise ValueError(f"invalid counts on {path} row {row_number}")
            if not math.isfinite(row_r2) or (row_r2 < 0.0 and row_r2 != -1.0) or row_r2 > 1.0:
                raise ValueError(f"invalid r-squared on {path} row {row_number}")
            if not math.isfinite(row_accuracy) or not 0.0 <= row_accuracy <= 1.0:
                raise ValueError(f"invalid accuracy on {path} row {row_number}")
            macs.append(mac)
            n_loci.append(loci)
            n_predictions.append(predictions)
            r2.append(row_r2)
            accuracy.append(row_accuracy)
    if not macs:
        raise ValueError(f"imputation evaluation contains no MAC bins: {path}")
    return macs, n_loci, n_predictions, r2, accuracy


def make_figure(
    macs: list[int],
    n_loci: list[int],
    n_predictions: list[int],
    r2: list[float],
    accuracy: list[float],
) -> go.Figure:
    custom = list(zip(n_loci, n_predictions))
    hover = (
        "reference MAC=%{x}<br>%{fullData.name}=%{y:.4f}"
        "<br>retained loci=%{customdata[0]:,}<br>target alleles=%{customdata[1]:,}<extra></extra>"
    )
    figure = go.Figure()
    figure.add_trace(go.Scatter(
        x=macs,
        y=[None if value < 0.0 else value for value in r2],
        customdata=custom,
        mode="lines",
        name="Pooled r²",
        line={"color": "#3366cc", "width": 2},
        hovertemplate=hover,
    ))
    figure.add_trace(go.Scatter(
        x=macs,
        y=accuracy,
        customdata=custom,
        mode="lines",
        name="Pooled accuracy",
        line={"color": "#ef8a35", "width": 2},
        hovertemplate=hover,
    ))
    total_predictions = sum(n_predictions)
    overall_accuracy = sum(n * value for n, value in zip(n_predictions, accuracy)) / total_predictions
    figure.update_layout(
        title="Imputation performance by reference minor-allele count",
        xaxis_title="Reference minor-allele count (haplotypes)",
        yaxis_title="Metric pooled across retained target alleles",
        yaxis_range=[0, 1],
        template="plotly_white",
        hovermode="x unified",
        legend={"orientation": "h", "y": 1.08, "x": 1, "xanchor": "right"},
        margin={"l": 70, "r": 30, "t": 90, "b": 65},
    )
    figure.add_annotation(
        x=0,
        y=1.08,
        xref="paper",
        yref="paper",
        showarrow=False,
        text=(
            f"retained loci={sum(n_loci):,}; target alleles={total_predictions:,}; "
            f"overall accuracy={overall_accuracy:.4f}"
        ),
    )
    return figure


def main() -> None:
    args = parse_args()
    values = read_evaluation(args.evaluation)
    figure = make_figure(*values)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(
        args.output,
        include_plotlyjs=ensure_plotly_asset(args.output),
        config={"responsive": True, "displaylogo": False},
    )
    print(
        f"wrote {args.output}: mac_bins={len(values[0])} retained_loci={sum(values[1])} "
        f"target_alleles={sum(values[2])}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
