#!/usr/bin/env python3

"""Write clean Plotly plots for imputation, blocks, and groups, plus IoU numbers.

Input run manifest is TSV with required columns ``panel``, ``k`` and
``fit_json``. Optional columns are ``positions``, ``init``, ``step1``,
``step3`` and ``converged``. Assignment paths are DFRA files. Position input
may be the comma-separated variant_pos.txt format or a VCF queried by
bcftools. Each fit JSON must contain ``block_snp_offsets`` and, when available,
``block_emission_counts``. ``--r2 LABEL=TSV`` accepts the standard imputation
evaluation TSV (mac, n_loci, n_predictions, r2). The output is a plain HTML
index linking three plots and a numerical partition-change summary.
"""

import argparse
import array
import csv
import json
import math
import struct
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from scipy.stats import gaussian_kde
from plotly.subplots import make_subplots
from plotly.colors import qualitative

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from plotly_html import ensure_plotly_asset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="?", help="optional TSV run manifest")
    parser.add_argument("--r2", action="append", default=[], metavar="LABEL=TSV",
                        help="imputation r-squared TSV; repeat for Beagle and DFCP runs")
    parser.add_argument("--genetic-map", type=Path,
                        default=Path("data/1000g_phase3_v5b/plink.chr20.GRCh37.map"),
                        help="PLINK or bp/cM map for genetic spans")
    parser.add_argument("--existing-root", type=Path,
                        default=Path("data/1000g_phase3_v5b/windows_64"))
    parser.add_argument("--maf1-root", type=Path,
                        default=Path("data/1000g_phase3_v5b/windows_8_maf1"))
    parser.add_argument("--output-dir", type=Path, default=Path("output/maf1_pbwt_study"))
    parser.add_argument("--output", type=Path, default=Path("output/maf1_pbwt_study/report.html"))
    return parser.parse_args()


def named_path(value: str, kind: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"{kind} must have the form LABEL=PATH: {value}")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise ValueError(f"invalid {kind}: {value}")
    return label, Path(raw_path)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    required = {"panel", "k", "fit_json"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"manifest must contain {', '.join(sorted(required))}: {path}")
    for number, row in enumerate(rows, start=2):
        if not row["panel"] or not row["fit_json"]:
            raise ValueError(f"missing panel or fit_json on manifest row {number}")
        try:
            row["k"] = str(int(row["k"]))
        except ValueError as error:
            raise ValueError(f"invalid K on manifest row {number}") from error
        for field in ("fit_json", "positions", "init", "step1", "step3", "converged",
                      "converged_fit"):
            if row.get(field):
                candidate = Path(row[field])
                row[field] = str(candidate if candidate.is_absolute() else path.parent / candidate)
    return rows


def read_positions(path: Path) -> list[int]:
    if path.suffix in {".vcf", ".gz", ".bgz"}:
        process = subprocess.run(
            ["bcftools", "query", "-f", "%POS\\n", str(path)],
            check=True, capture_output=True, text=True,
        )
        positions = [int(value) for value in process.stdout.split()]
    else:
        fields = path.read_text().replace(",", " ").split()
        try:
            declared = int(fields[0])
            positions = [int(value) for value in fields[1:]]
        except (IndexError, ValueError) as error:
            raise ValueError(f"invalid variant position file: {path}") from error
        if declared != len(positions):
            raise ValueError(f"position count does not match contents: {path}")
    if not positions or any(b < a for a, b in zip(positions, positions[1:])):
        raise ValueError(f"positions must be nonempty and nondecreasing: {path}")
    return positions


def read_map(path: Path) -> tuple[list[int], list[float]]:
    positions: list[int] = []
    cm: list[float] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            if len(fields) >= 4:  # PLINK: chromosome, marker, cM, bp
                bp, genetic = int(fields[3]), float(fields[2])
            elif len(fields) == 3:  # bp, cM/Mb, cM
                bp, genetic = int(fields[0]), float(fields[2])
            elif len(fields) == 2:  # bp, cM
                bp, genetic = int(fields[0]), float(fields[1])
            else:
                continue
        except ValueError:
            continue
        if positions and bp <= positions[-1]:
            raise ValueError(f"genetic map bp positions must increase: {path}")
        if cm and genetic < cm[-1]:
            raise ValueError(f"genetic map cM positions must not decrease: {path}")
        positions.append(bp)
        cm.append(genetic)
    if len(positions) < 2:
        raise ValueError(f"genetic map has fewer than two usable rows: {path}")
    return positions, cm


def interpolate_cm(position: int, map_bp: list[int], map_cm: list[float]) -> float:
    import bisect
    right = min(max(bisect.bisect_right(map_bp, position), 1), len(map_bp) - 1)
    left = right - 1
    fraction = (position - map_bp[left]) / (map_bp[right] - map_bp[left])
    return map_cm[left] + fraction * (map_cm[right] - map_cm[left])


def read_fit(path: Path) -> dict:
    value = json.loads(path.read_text())
    # Shared fits keep these fields at the root; target-fit JSON may wrap runs.
    if "block_snp_offsets" not in value and value.get("target_runs"):
        value = value["target_runs"][0]
    return value


def add_standard_runs(rows: list[dict[str, str]], root: Path, panel: str,
                      run_names: list[str], k_by_run: dict[str, int],
                      stages: dict[str, str] | None = None) -> None:
    """Expand the repository's per-window impute JSON naming convention."""
    for name in run_names:
        for fit_path in sorted(root.glob(f"window_*/impute/{name}.json")):
            window = fit_path.parents[1]
            row = {
                "panel": panel,
                "k": str(k_by_run[name]),
                "fit_json": str(fit_path),
                "positions": str(window / "ref.vcf.gz"),
                "init": "", "step1": "", "step3": "", "converged": "",
            }
            if stages and name in stages:
                stage_prefix = stages[name]
                suffixes = {"init": "step0", "step1": "step1", "step3": "step3",
                            "converged": "converged"}
                for field, suffix in suffixes.items():
                    assign_path = window / "impute" / f"{stage_prefix}_{suffix}.r_assign.bin"
                    if assign_path.is_file():
                        row[field] = str(assign_path)
                converged_fit = window / "impute" / f"{stage_prefix}_converged.json"
                if converged_fit.is_file():
                    row["converged_fit"] = str(converged_fit)
                    row["converged_cap"] = "30"
            rows.append(row)


def default_manifest(args: argparse.Namespace) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    existing_names = [f"joint_k{k}_step1" for k in (256, 512, 1024, 2048)]
    maf_names = [f"maf1_joint_k{k}_step1" for k in (128, 256, 512, 1024, 2048)]
    add_standard_runs(rows, args.existing_root, "existing", existing_names,
                      {name: int(name.split("k")[1].split("_")[0]) for name in existing_names},
                      {"joint_k1024_step1": "joint_k1024"})
    add_standard_runs(rows, args.maf1_root, "MAF ≥1%", maf_names,
                      {name: int(name.split("k")[1].split("_")[0]) for name in maf_names},
                      {"maf1_joint_k1024_step1": "maf1_joint_k1024"})

    if not args.r2:
        output = args.output_dir
        for k in (256, 512, 1024, 2048):
            path = Path(f"output/joint_k{k}_step1.tsv")
            if path.is_file():
                args.r2.append(f"existing K{k}={path}")
        for k in (128, 256, 512, 1024, 2048):
            path = output / f"maf1_joint_k{k}_step1.tsv"
            if path.is_file() and k != 1024:
                args.r2.append(f"MAF ≥1% K{k}={path}")
        for stage in ("step0", "step1", "step3", "converged"):
            path = output / f"maf1_joint_k1024_{stage}.tsv"
            if path.is_file():
                args.r2.append(f"MAF K1024 {stage}={path}")
        beagle_paths = (
            ("Beagle existing", Path("data/1000g_phase3_v5b/windows/eval_impute_beagle.tsv")),
            ("Beagle MAF ≥1%", Path("data/1000g_phase3_v5b/maf1/eval_impute_beagle.tsv")),
        )
        for label, path in beagle_paths:
            if path.is_file():
                args.r2.insert(0, f"{label}={path}")
    return rows


def read_dfra(path: Path) -> tuple[int, int, array.array]:
    with path.open("rb") as stream:
        header = stream.read(12)
        if len(header) != 12 or header[:4] != b"DFRA":
            raise ValueError(f"invalid DFRA assignment file: {path}")
        n, loci = struct.unpack_from("<II", header, 4)
        if n == 0 or loci == 0 or path.stat().st_size != 12 + n * loci * 4:
            raise ValueError(f"DFRA dimensions do not match file size: {path}")
        ids = array.array("I")
        ids.fromfile(stream, n * loci)
    if ids.itemsize != 4:
        raise RuntimeError("platform unsigned int is not 32 bits")
    if sys.byteorder != "little":
        ids.byteswap()
    return n, loci, ids


def read_seq_n(path: Path) -> int:
    with path.open("rb") as stream:
        header = stream.read(12)
    if len(header) != 12 or header[:4] != b"DFCP":
        raise ValueError(f"invalid DFCP reference file: {path}")
    n, loci = struct.unpack_from("<II", header, 4)
    if n == 0 or loci == 0:
        raise ValueError(f"empty DFCP reference file: {path}")
    return n


def choose_two(value: int) -> int:
    return value * (value - 1) // 2


def partition_iou(left: list[int], right: list[int]) -> float:
    left_counts = Counter(left)
    right_counts = Counter(right)
    joint_counts = Counter(zip(left, right))
    intersection = sum(choose_two(count) for count in joint_counts.values())
    union = (sum(choose_two(count) for count in left_counts.values())
             + sum(choose_two(count) for count in right_counts.values()) - intersection)
    return 1.0 if union == 0 else intersection / union


def compare_assignments(init_path: Path, stage_path: Path, offsets: list[int]) -> list[float]:
    n_init, l_init, init = read_dfra(init_path)
    n_stage, l_stage, stage = read_dfra(stage_path)
    if (n_init, l_init) != (n_stage, l_stage):
        raise ValueError(f"assignment dimensions differ: {init_path} and {stage_path}")
    if offsets[0] != 0 or offsets[-1] != l_init or any(b <= a for a, b in zip(offsets, offsets[1:])):
        raise ValueError(f"block offsets do not partition the assignment loci: {init_path}")
    scores = []
    for start, end in zip(offsets, offsets[1:]):
        left = [init[i * l_init + start] for i in range(n_init)]
        right = [stage[i * l_stage + start] for i in range(n_stage)]
        scores.append(partition_iou(left, right))
    return scores


def summarize_run(row: dict[str, str], map_data: tuple[list[int], list[float]] | None) -> tuple[dict, dict]:
    fit = read_fit(Path(row["fit_json"]))
    offsets = fit.get("block_snp_offsets")
    counts = fit.get("block_emission_counts")
    block = {"panel": row["panel"], "k": int(row["k"]), "label": f"{row['panel']} K={row['k']}"}
    if offsets:
        block["summary_level"] = "per-block"
        if any(b <= a for a, b in zip(offsets, offsets[1:])):
            raise ValueError(f"block offsets are not increasing: {row['fit_json']}")
        snp_lengths = [b - a for a, b in zip(offsets, offsets[1:])]
        block["snp_lengths"] = snp_lengths
        if counts:
            if len(counts) != len(snp_lengths):
                raise ValueError(f"emission counts do not match offsets: {row['fit_json']}")
            if any(count <= 0 for count in counts):
                raise ValueError(f"emission counts must be positive: {row['fit_json']}")
            block["emission_counts"] = counts
            n_sequences = int(fit.get("n_ref_selected", fit.get("n_ref", 0)))
            if not n_sequences and row.get("init"):
                n_sequences = read_dfra(Path(row["init"]))[0]
            if not n_sequences and fit.get("ref_file"):
                n_sequences = read_seq_n(Path(fit["ref_file"]))
            if n_sequences:
                block["mean_group_sizes"] = [n_sequences * len(counts) / sum(counts)]
                block["group_size_level"] = "per-window mean across emissions"
        if row.get("positions"):
            positions = read_positions(Path(row["positions"]))
            if offsets[-1] != len(positions):
                raise ValueError(f"position count differs from block offsets: {row['positions']}")
            block["bp_lengths"] = [positions[end - 1] - positions[start]
                                   for start, end in zip(offsets, offsets[1:])]
            if map_data:
                map_bp, map_cm = map_data
                block["cm_lengths"] = [
                    interpolate_cm(positions[end - 1], map_bp, map_cm)
                    - interpolate_cm(positions[start], map_bp, map_cm)
                    for start, end in zip(offsets, offsets[1:])
                ]
    else:
        # Keep old runs usable as single summary points.
        block["summary_level"] = "per-window mean (legacy JSON)"
        for field, key in (("snp_lengths", "mean_block_snps"),
                           ("emission_counts", "mean_block_emissions")):
            if key in fit:
                block[field] = [float(fit[key])]
        if "mean_block_emissions" in fit and fit.get("ref_file"):
            n_sequences = read_seq_n(Path(fit["ref_file"]))
            block["mean_group_sizes"] = [n_sequences / float(fit["mean_block_emissions"])]
            block["group_size_level"] = "per-window mean across emissions"

    stages = {}
    if row.get("converged_fit"):
        converged_fit = read_fit(Path(row["converged_fit"]))
        converged_steps = len(converged_fit.get("train_log", []))
        converged_cap = int(row.get("converged_cap", 30))
        block["converged_steps"] = converged_steps
        block["converged_cap"] = converged_cap
        block["converged_cap_hit"] = converged_steps >= converged_cap
    init_path = row.get("init")
    if init_path and offsets:
        init_assignment = Path(init_path)
        init_json = init_assignment.with_name(
            init_assignment.name.removesuffix(".r_assign.bin") + ".json"
        )
        if not init_json.is_file():
            raise ValueError(f"initial fit JSON is missing for assignment: {init_assignment}")
        if read_fit(init_json).get("block_snp_offsets") != offsets:
            raise ValueError(f"block offsets differ between {row['fit_json']} and {init_json}")
        for stage in ("step1", "step3", "converged"):
            if row.get(stage):
                stage_assignment = Path(row[stage])
                stage_json = stage_assignment.with_name(
                    stage_assignment.name.removesuffix(".r_assign.bin") + ".json"
                )
                if not stage_json.is_file():
                    raise ValueError(f"stage fit JSON is missing for assignment: {stage_assignment}")
                stage_offsets = read_fit(stage_json).get("block_snp_offsets")
                if stage_offsets != offsets:
                    raise ValueError(f"block offsets differ between {row['fit_json']} and {stage_json}")
                stages[stage] = compare_assignments(
                    Path(init_path), Path(row[stage]), offsets
                )
    return block, stages


def read_r2(value: str) -> tuple[str, list[tuple[int, int, int, float]]]:
    label, path = named_path(value, "--r2")
    rows: list[tuple[int, int, int, float]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"mac", "n_loci", "n_predictions", "r2"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"unexpected imputation TSV header: {path}")
        for number, row in enumerate(reader, start=2):
            mac = int(row["mac"])
            n_loci = int(row["n_loci"])
            n_predictions = int(row["n_predictions"])
            r2 = float(row["r2"])
            if (mac < 0 or n_loci <= 0 or n_predictions <= 0
                    or not math.isfinite(r2) or r2 > 1 or r2 < -1):
                raise ValueError(f"invalid evaluation values on {path} row {number}")
            rows.append((mac, n_loci, n_predictions, r2))
    return label, rows


def plot_imputation(runs: list[tuple[str, list[tuple[int, int, int, float]]]]) -> go.Figure:
    figure = go.Figure()
    filtered_ks = (128, 256, 512, 1024, 2048)
    colors = {k: qualitative.Plotly[i] for i, k in enumerate(filtered_ks)}
    for label, values in runs:
        is_filtered_k = label.startswith("MAF ≥1% K")
        if label != "Beagle MAF ≥1%" and not is_filtered_k:
            continue
        line = {"width": 2, "color": "#111111" if label == "Beagle MAF ≥1%" else None}
        if is_filtered_k:
            line["color"] = colors[int(label.rsplit("K", 1)[1])]
        figure.add_trace(go.Scatter(
            x=[mac for mac, _, _, _ in values],
            y=[None if r2 < 0 else r2 for _, _, _, r2 in values],
            customdata=[loci for _, loci, _, _ in values],
            mode="lines", name=label, line=line,
            hovertemplate="MAC=%{x}<br>r²=%{y:.4f}<br>loci=%{customdata:,}<extra>%{fullData.name}</extra>",
        ))
    figure.update_layout(
        title="Imputation performance by reference minor-allele count",
        xaxis={"title": "Reference minor-allele count (haplotypes)", "type": "log"},
        yaxis={"title": "Imputation r²", "range": [0, 1]},
        hovermode="x unified", height=760, template="plotly_white",
        legend={"x": 1.02, "y": 1, "xanchor": "left"},
        margin={"l": 75, "r": 280, "t": 90, "b": 65},
    )
    return figure


def density_curve(values: list[float]) -> tuple[np.ndarray, np.ndarray]:
    """KDE in log(1 + length), so long-tailed block lengths remain readable."""
    transformed = np.log1p(np.asarray(values, dtype=float))
    upper = max(float(np.quantile(transformed, 0.995)), 1e-6)
    grid = np.linspace(0, upper, 300)
    if len(transformed) < 2 or np.allclose(transformed, transformed[0]):
        center = float(transformed[0])
        return grid, np.exp(-0.5 * ((grid - center) / max(upper / 100, 1e-6)) ** 2)
    return grid, gaussian_kde(transformed)(grid)


def plot_lengths(blocks: list[dict]) -> go.Figure:
    figure = go.Figure()
    units = (("snp_lengths", "SNPs"), ("bp_lengths", "bp"), ("cm_lengths", "cM"))
    filtered_ks = (128, 256, 512, 1024, 2048)
    visible_by_unit: dict[str, list[bool]] = {label: [] for _, label in units}
    for key, unit in units:
        for i, k in enumerate(filtered_ks):
            values = [value for block in blocks if block["panel"] == "MAF ≥1%" and block["k"] == k
                      for value in block.get(key, [])]
            if not values:
                continue
            grid, density = density_curve(values)
            figure.add_trace(go.Scatter(
                x=grid, y=density, mode="lines", name=f"K={k}",
                line={"width": 3, "color": qualitative.Plotly[i]},
                visible=(unit == "SNPs"),
                customdata=np.expm1(grid),
                hovertemplate=f"K={k}<br>length=%{{customdata:.3g}} {unit}<br>density=%{{y:.3g}}<extra></extra>",
            ))
            for label in visible_by_unit:
                visible_by_unit[label].append(label == unit)
    ticks = {
        "SNPs": [0, 1, 2, 5, 10, 25, 50, 100, 250, 500, 1000],
        "bp": [0, 100, 1000, 10000, 100000, 1000000],
        "cM": [0, 0.001, 0.01, 0.1, 1, 10],
    }
    buttons = []
    for unit in visible_by_unit:
        buttons.append(dict(
            label=unit, method="update",
            args=[{"visible": visible_by_unit[unit]}, {
                "xaxis.tickvals": np.log1p(ticks[unit]).tolist(),
                "xaxis.ticktext": [str(tick) for tick in ticks[unit]],
                "xaxis.title.text": f"Block length ({unit}; logarithmic spacing)",
            }],
        ))
    figure.update_layout(
        title="Greedy PBWT block-length density (MAF ≥1% panel)",
        xaxis={"title": "Block length (SNPs; logarithmic spacing)",
               "tickvals": np.log1p(ticks["SNPs"]).tolist(),
               "ticktext": [str(tick) for tick in ticks["SNPs"]]},
        yaxis={"title": "Density in log(1 + length)", "rangemode": "tozero"},
        updatemenus=[{"buttons": buttons, "direction": "right", "type": "buttons",
                      "x": 1.02, "y": 0.45, "xanchor": "left", "showactive": True}],
        template="plotly_white", height=760,
        margin={"l": 75, "r": 280, "t": 90, "b": 65},
        legend={"x": 1.02, "y": 1, "xanchor": "left"},
    )
    return figure


def plot_groups(blocks: list[dict]) -> go.Figure:
    figure = make_subplots(rows=2, cols=1, vertical_spacing=0.15,
                           subplot_titles=("Mean exact emission patterns per block",
                                           "Mean reference haplotypes per pattern"))
    ks = (128, 256, 512, 1024, 2048)
    for panel, color, dash in (("existing", "#6B7280", "dash"),
                               ("MAF ≥1%", "#2563EB", "solid")):
        for row, key, hover in ((1, "emission_counts", "patterns/block"),
                                (2, "mean_group_sizes", "haplotypes/pattern")):
            xs = []
            ys = []
            for k in ks:
                values = [value for block in blocks if block["panel"] == panel and block["k"] == k
                          for value in block.get(key, [])]
                if values:
                    xs.append(k)
                    ys.append(float(np.mean(values)))
            figure.add_trace(go.Scatter(
                x=xs, y=ys, mode="lines+markers", name=panel, legendgroup=panel,
                showlegend=row == 1,
                line={"color": color, "width": 2, "dash": dash}, marker={"size": 8},
                hovertemplate=f"K=%{{x}}<br>{hover}=%{{y:.2f}}<extra>%{{fullData.name}}</extra>",
            ), row=row, col=1)
    figure.update_layout(title="PBWT pattern groups", height=850, template="plotly_white",
                         margin={"l": 75, "r": 280, "t": 90, "b": 65},
                         legend={"x": 1.02, "y": 1, "xanchor": "left", "groupclick": "togglegroup"})
    figure.update_yaxes(title_text="Patterns per block", row=1, col=1)
    figure.update_yaxes(title_text="Haplotypes per pattern", row=2, col=1)
    figure.update_xaxes(title_text="Maximum patterns K", row=1, col=1)
    figure.update_xaxes(title_text="Maximum patterns K", row=2, col=1)
    return figure


def partition_rows(stage_runs: list[tuple[dict, dict]]) -> list[tuple[str, int, float, float]]:
    values_by_stage = {stage: [] for stage in ("step1", "step3", "converged")}
    for block, stages in stage_runs:
        if block["panel"] == "MAF ≥1%" and block["k"] == 1024:
            for stage, values in stages.items():
                values_by_stage[stage].extend(values)
    labels = {"step1": "1 step", "step3": "3 steps", "converged": "Converged"}
    return [(labels[stage], len(values), float(np.mean(values)), float(np.median(values)))
            for stage, values in values_by_stage.items() if values]


def write_partition_table(path: Path, rows: list[tuple[str, int, float, float]]) -> None:
    body = "".join(
        f"<tr><td>{stage}</td><td>{count:,}</td><td>{mean:.4f}</td><td>{median:.4f}</td></tr>"
        for stage, count, mean, median in rows
    )
    path.write_text(
        '<!doctype html><html><head><meta charset="utf-8"><title>Partition change</title>'
        '<style>body{font:16px system-ui,sans-serif;color:#222;background:white;max-width:900px;'
        'margin:48px auto;line-height:1.5}table{border-collapse:collapse;width:100%}'
        'th,td{text-align:left;padding:10px 14px;border-bottom:1px solid #ddd}'
        'th{font-weight:600}a{color:#2459a6}</style></head><body>'
        '<h1>Partition change from initialization</h1>'
        '<p>K=1024, eight MAF ≥1% windows. Pairwise co-clustering IoU is label-invariant; '
        '1 means the same partition.</p><table><thead><tr><th>Stage</th><th>Blocks</th>'
        '<th>Mean IoU</th><th>Median IoU</th></tr></thead><tbody>' + body + '</tbody></table>'
        '<p>Convergence took 7–11 steps; none of the eight windows reached the 30-step cap.</p>'
        '</body></html>'
    )


def write_index(path: Path, rows: list[tuple[str, int, float, float]]) -> None:
    links = (("Imputation accuracy", "imputation.html"),
             ("Block-length density", "block_lengths.html"),
             ("Emission and group sizes", "groups.html"),
             ("Partition change numbers", "partition_iou.html"))
    items = "".join(f'<li><a href="{filename}">{title}</a></li>' for title, filename in links)
    body = "".join(f"<tr><td>{stage}</td><td>{mean:.4f}</td><td>{median:.4f}</td></tr>"
                   for stage, _, mean, median in rows)
    path.write_text(
        '<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" '
        'content="width=device-width, initial-scale=1"><title>PBWT study</title><style>'
        'body{font:16px system-ui,sans-serif;color:#222;background:white;max-width:900px;'
        'margin:48px auto;line-height:1.5}a{color:#2459a6}'
        'table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:8px 12px;'
        'border-bottom:1px solid #ddd}</style></head><body><h1>PBWT study</h1><ul>'
        + items + '</ul><h2>Partition change from initialization (K=1024)</h2>'
        '<table><tr><th>Stage</th><th>Mean IoU</th><th>Median IoU</th></tr>'
        + body + '</table></body></html>'
    )


def main() -> None:
    args = parse_args()
    if args.runs:
        rows = read_manifest(args.runs)
    else:
        rows = default_manifest(args)
    map_data = read_map(args.genetic_map) if args.genetic_map else None
    blocks = []
    stage_runs = []
    for row in rows:
        block, stages = summarize_run(row, map_data)
        blocks.append(block)
        if stages:
            stage_runs.append((block, stages))
    r2_runs = [read_r2(value) for value in args.r2]
    maf_beagle = next((values for label, values in r2_runs if label == "Beagle MAF ≥1%"), None)
    maf_k128 = next((values for label, values in r2_runs if label == "MAF ≥1% K128"), None)
    if maf_beagle is not None and maf_k128 is not None:
        beagle_totals = (sum(row[1] for row in maf_beagle), sum(row[2] for row in maf_beagle))
        dfcp_totals = (sum(row[1] for row in maf_k128), sum(row[2] for row in maf_k128))
        if beagle_totals != dfcp_totals:
            raise ValueError(f"filtered Beagle and K128 evaluation totals differ: {beagle_totals} != {dfcp_totals}")
        print(f"verified filtered Beagle/K128 evaluation totals: loci={beagle_totals[0]} predictions={beagle_totals[1]}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figures = {
        "imputation.html": plot_imputation(r2_runs),
        "block_lengths.html": plot_lengths(blocks),
        "groups.html": plot_groups(blocks),
    }
    for filename, figure in figures.items():
        path = args.output.parent / filename
        figure.write_html(path, include_plotlyjs=ensure_plotly_asset(path),
                          config={"responsive": True, "displaylogo": False})
    rows = partition_rows(stage_runs)
    write_partition_table(args.output.parent / "partition_iou.html", rows)
    write_index(args.output, rows)
    print(f"wrote {args.output}, {len(figures)} plots, and partition numbers")


if __name__ == "__main__":
    main()
