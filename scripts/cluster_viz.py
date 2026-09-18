#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.colors import qualitative
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter1d

from plotly_html import ensure_plotly_asset


CLADE_FIELDS = ["locus", "position", "r_id", "size", "weight", "iou", "time"]
TRACT_FIELDS = [
    "start_locus", "end_locus", "start_position", "end_position", "size", "n_loci", "span_bp",
]


def read_values(
    path: Path, fields: list[str], value_field: str, weighted: bool
) -> tuple[np.ndarray, np.ndarray]:
    values = []
    weights = []
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if rows.fieldnames != fields:
            raise ValueError(f"unexpected cluster evaluation header: {path}")
        for row in rows:
            weight = float(row["weight"]) if weighted else 1.0
            if not weighted or weight > 0.0:
                values.append(float(row[value_field]))
                weights.append(weight)
    if not values:
        raise ValueError(f"cluster evaluation has no values: {path}")
    return np.asarray(values), np.asarray(weights)


def densities(
    series: list[tuple[np.ndarray, np.ndarray]], log_scale: bool = False
) -> tuple[np.ndarray, list[np.ndarray]]:
    transformed = [(np.log1p(values) if log_scale else values, weights) for values, weights in series]
    maximum = max(1.0, max(values.max() for values, _ in transformed))
    edges = np.linspace(0.0, maximum, 257)
    centers = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]
    result = []
    for values, weights in transformed:
        density = np.histogram(values, bins=edges, weights=weights)[0] / (weights.sum() * width)
        result.append(gaussian_filter1d(density, sigma=2.0))
    if log_scale:
        centers = np.expm1(centers)
    return centers, result


def display_name(name: str) -> str:
    if name == "beagle4":
        return "Beagle 4 DAG"
    if match := re.fullmatch(r"pbwt(\d+)(?:_(init|step\d+|converged))?", name):
        stage = f", {match.group(2)}" if match.group(2) else ""
        return f"PBWT r={match.group(1)}{stage}"
    if match := re.fullmatch(r"greedy(\d+)(?:_(init|step\d+|converged))?", name):
        stage = f", {match.group(2)}" if match.group(2) else ""
        return f"Greedy K={match.group(1)}{stage}"
    return name.replace("_", " ")


def make_figure(
    evaluations: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
) -> go.Figure:
    time_series = [(time, time_weight) for _, time, time_weight, _, _, _, _ in evaluations]
    time_x, time_densities = densities(time_series, log_scale=True)
    tract_bp_series = [(tract, weight) for _, _, _, tract, weight, _, _ in evaluations]
    tract_bp_x, tract_bp_densities = densities(tract_bp_series, log_scale=True)
    tract_loci_series = [(tract, weight) for _, _, _, _, _, tract, weight in evaluations]
    tract_loci_x, tract_loci_densities = densities(tract_loci_series, log_scale=True)
    figure = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Best-matching true-clade time", "Exact-membership R-cluster tract length"),
        vertical_spacing=0.13,
    )
    for i, (label, _, _, _, _, _, _) in enumerate(evaluations):
        color = qualitative.Plotly[i % len(qualitative.Plotly)]
        figure.add_trace(go.Scatter(
            x=time_x,
            y=time_densities[i],
            mode="lines",
            name=label,
            legendgroup=label,
            line={"color": color, "width": 2},
            hovertemplate="time=%{x:.3g} generations<br>density=%{y:.3g}<extra>%{fullData.name}</extra>",
        ), row=1, col=1)
        figure.add_trace(go.Scatter(
            x=tract_bp_x,
            y=tract_bp_densities[i],
            mode="lines",
            name=label,
            legendgroup=label,
            line={"color": color, "width": 2},
            showlegend=False,
            hovertemplate="span=%{x:.3g} bp<br>density=%{y:.3g}<extra>%{fullData.name}</extra>",
        ), row=2, col=1)
    figure.update_layout(
        title="DFCP cluster evaluation densities",
        template="plotly_white",
        legend={"y": 1, "x": 1.02, "xanchor": "left", "groupclick": "togglegroup"},
        margin={"l": 75, "r": 280, "t": 90, "b": 65},
        updatemenus=[{
            "type": "buttons",
            "direction": "right",
            "x": 1.02,
            "y": 0.45,
            "xanchor": "left",
            "yanchor": "middle",
            "buttons": [
                {
                    "label": "bp",
                    "method": "update",
                    "args": [
                        {
                            "x": [tract_bp_x] * len(evaluations),
                            "y": tract_bp_densities,
                            "hovertemplate": [
                                "span=%{x:.3g} bp<br>density=%{y:.3g}<extra>%{fullData.name}</extra>"
                            ] * len(evaluations),
                        },
                        {
                            "xaxis2.title.text": "Physical tract span (bp; log scale)",
                            "yaxis2.title.text": "Tract density in log-span",
                        },
                        list(range(1, 2 * len(evaluations), 2)),
                    ],
                },
                {
                    "label": "loci",
                    "method": "update",
                    "args": [
                        {
                            "x": [tract_loci_x] * len(evaluations),
                            "y": tract_loci_densities,
                            "hovertemplate": [
                                "length=%{x:.3g} loci<br>density=%{y:.3g}<extra>%{fullData.name}</extra>"
                            ] * len(evaluations),
                        },
                        {
                            "xaxis2.title.text": "Exact tract length (loci; log scale)",
                            "yaxis2.title.text": "Tract density in log-length",
                        },
                        list(range(1, 2 * len(evaluations), 2)),
                    ],
                },
            ],
        }],
    )
    figure.update_xaxes(
        title_text="Clade time (generations before present; log scale)", type="log", row=1, col=1
    )
    figure.update_xaxes(title_text="Physical tract span (bp; log scale)", type="log", row=2, col=1)
    figure.update_yaxes(title_text="Clade-weighted density in log-time", row=1, col=1)
    figure.update_yaxes(title_text="Tract density in log-span", row=2, col=1)
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare DFCP cluster evaluation densities.")
    parser.add_argument("evaluations", type=Path, nargs="+", help="eval_clusters JSON output files")
    parser.add_argument("--output", type=Path, default=Path("clusters.html"))
    args = parser.parse_args()

    evaluations = []
    for path in args.evaluations:
        suffix = ".eval.json"
        if not path.name.endswith(suffix):
            raise ValueError(f"evaluation filename must end in {suffix}: {path}")
        base = path.with_name(path.name[:-len(suffix)])
        clade_times = read_values(
            base.with_name(base.name + ".clade_times.tsv"), CLADE_FIELDS, "time", True
        )
        tract_bp = read_values(
            base.with_name(base.name + ".cluster_tracts.tsv"), TRACT_FIELDS, "span_bp", False
        )
        tract_loci = read_values(
            base.with_name(base.name + ".cluster_tracts.tsv"), TRACT_FIELDS, "n_loci", False
        )
        evaluations.append((display_name(base.name), *clade_times, *tract_bp, *tract_loci))
    figure = make_figure(evaluations)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(
        args.output,
        include_plotlyjs=ensure_plotly_asset(args.output),
        config={"responsive": True, "displaylogo": False},
    )
    print(f"wrote {args.output}: evaluations={len(evaluations)}")


if __name__ == "__main__":
    main()
