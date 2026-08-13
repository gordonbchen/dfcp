"""Build an interactive report for one DFCP JSON result."""

import argparse
import html
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from plotly_html import ensure_plotly_asset


COLORS = {
    "dfcp": "#5B5BD6",
    "emission": "#E8794F",
    "green": "#28A17B",
    "gold": "#D5A52E",
    "ink": "#243447",
}

DISPLAY_NAMES = {
    "viterbi_impute_acc": "Viterbi",
    "fwd_bkwd_impute_acc": "Forward-backward",
    "mode_impute_acc": "Locus mode",
    "mean_iou": "DFCP adjacent-locus IoU",
    "mean_emission_iou": "Emission adjacent-locus IoU",
    "clade_iou": "DFCP clade IoU",
    "emission_clade_iou": "Emission clade IoU",
    "mean_excess_parsimony": "DFCP excess parsimony",
    "mean_emission_excess_parsimony": "Emission excess parsimony",
    "mean_clusters": "Mean clusters / locus",
    "cluster_purity": "Cluster purity",
    "n_val_seqs": "Held-out sequences",
    "n_masked_ls": "Masked loci",
    "clade_beta": "Clade weight beta",
    "t_init": "Initialization",
    "t_viterbi_impute": "Viterbi imputation",
    "t_fwd_bkwd_impute": "Forward-backward imputation",
    "t_parsimony": "Parsimony",
    "t_clade_iou": "Clade IoU",
    "t_iou": "Adjacent-locus IoU",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an interactive Plotly report from one DFCP JSON result."
    )
    parser.add_argument(
        "json", type=Path, nargs="?",
        help="result JSON (reads stdin when omitted)",
    )
    parser.add_argument(
        "--output", "--out_file", dest="output", type=Path,
        default=Path("viz.html"), help="HTML report path (default: viz.html)",
    )
    parser.add_argument(
        "--tree-dot", type=Path,
        help="optional DOT file produced by DFCP --tree_vis",
    )
    return parser.parse_args()


def load_result(path: Path | None) -> dict[str, Any]:
    if path is None:
        result = json.load(sys.stdin)
        source = "stdin"
    else:
        result = json.loads(path.read_text())
        source = str(path)
    if not isinstance(result, dict) or not isinstance(result.get("seq_file"), str):
        raise ValueError(f"{source} is not a DFCP result JSON object")
    return result


def display_name(name: str) -> str:
    return DISPLAY_NAMES.get(name, name.replace("_", " ").capitalize())


def format_value(value: int | float) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, int):
        return f"{value:,}"
    if abs(value) >= 10_000 or (value != 0 and abs(value) < 0.001):
        return f"{value:.4g}"
    return f"{value:.4f}".rstrip("0").rstrip(".")


def style_figure(fig: go.Figure, height: int = 390) -> go.Figure:
    fig.update_layout(
        height=height,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#FBFCFE",
        font={"family": "Inter, ui-sans-serif, system-ui, sans-serif", "color": COLORS["ink"]},
        margin={"l": 58, "r": 24, "t": 58, "b": 50},
        hoverlabel={"bgcolor": "white"},
        legend={"orientation": "h", "y": 1.13, "x": 0},
    )
    fig.update_xaxes(gridcolor="#E8EDF3", zeroline=False)
    fig.update_yaxes(gridcolor="#E8EDF3", zeroline=False)
    return fig


def elbo_figure(result: dict[str, Any]) -> go.Figure | None:
    train_log = result.get("train_log")
    if not isinstance(train_log, list) or not train_log:
        return None
    steps = list(range(1, len(train_log) + 1))
    fig = go.Figure(
        go.Scatter(
            x=steps, y=[step["elbo"] for step in train_log], mode="lines+markers",
            name="ELBO", showlegend=False, line={"color": COLORS["dfcp"], "width": 3},
            hovertemplate="iteration %{x}<br>ELBO %{y:.6g}<extra></extra>",
        )
    )
    fig.update_layout(title="ELBO over training")
    fig.update_xaxes(title_text="iteration", dtick=1)
    fig.update_yaxes(title_text="ELBO")
    return style_figure(fig, 360)


def iteration_time_figure(result: dict[str, Any]) -> go.Figure | None:
    train_log = result.get("train_log")
    if not isinstance(train_log, list) or not train_log:
        return None
    steps = list(range(1, len(train_log) + 1))
    fig = go.Figure()
    for key, color in (
        ("t_max", COLORS["dfcp"]),
        ("t_expect", COLORS["green"]),
        ("t_elbo", COLORS["gold"]),
        ("t_step", COLORS["ink"]),
    ):
        if all(key in step for step in train_log):
            fig.add_trace(
                go.Scatter(
                    x=steps, y=[step[key] for step in train_log], mode="lines+markers",
                    name=display_name(key), line={"color": color, "width": 2},
                    hovertemplate=f"iteration %{{x}}<br>{display_name(key)} %{{y:,}} ms<extra></extra>",
                )
            )
    fig.update_layout(title="Time per training iteration")
    fig.update_xaxes(title_text="iteration", dtick=1)
    fig.update_yaxes(title_text="milliseconds")
    return style_figure(fig, 360)


def parameter_figure(result: dict[str, Any]) -> go.Figure | None:
    params = result.get("params")
    if not isinstance(params, dict):
        return None
    series = [("mu_gamma", COLORS["emission"]), ("mu_d", COLORS["dfcp"])]
    available = [(key, color) for key, color in series if isinstance(params.get(key), list)]
    if not available:
        return None
    fig = make_subplots(
        rows=len(available), cols=1, shared_xaxes=True,
        subplot_titles=tuple(display_name(key) for key, _ in available),
        vertical_spacing=0.18,
    )
    for row, (key, color) in enumerate(available, start=1):
        values = params[key]
        fig.add_trace(
            go.Scattergl(
                x=list(range(len(values))), y=values, mode="lines",
                name=display_name(key), line={"color": color, "width": 2},
                hovertemplate=f"locus %{{x}}<br>{display_name(key)} %{{y:.5g}}<extra></extra>",
            ), row=row, col=1,
        )
        fig.update_yaxes(title_text=key.removeprefix("mu_"), row=row, col=1)
    fig.update_xaxes(title_text="locus", row=len(available), col=1)
    return style_figure(fig, 260 + 180 * len(available))


def imputation_figure(result: dict[str, Any]) -> go.Figure | None:
    specs = (
        ("viterbi_impute_acc", "t_viterbi_impute", COLORS["dfcp"]),
        ("fwd_bkwd_impute_acc", "t_fwd_bkwd_impute", COLORS["green"]),
        ("mode_impute_acc", None, COLORS["gold"]),
    )
    rows = [(acc, timing, color) for acc, timing, color in specs if acc in result]
    if not rows:
        return None
    labels = [display_name(acc) for acc, _, _ in rows]
    times = [
        f"{result[timing]:,} ms" if timing and timing in result else "not timed"
        for _, timing, _ in rows
    ]
    fig = go.Figure(go.Bar(
        x=labels,
        y=[result[acc] for acc, _, _ in rows],
        marker_color=[color for _, _, color in rows],
        text=[f"{result[acc]:.2%}" for acc, _, _ in rows],
        textposition="outside",
        customdata=times,
        hovertemplate=(
            "%{x}<br>accuracy %{y:.3%}<br>elapsed "
            "%{customdata}<extra></extra>"
        ),
    ))
    fig.update_layout(title="Held-out imputation accuracy")
    fig.update_yaxes(title_text="accuracy", range=[0, 1.08], tickformat=".0%")
    return style_figure(fig)


def metric_figure(result: dict[str, Any]) -> go.Figure | None:
    specs = (
        ("mean_iou", COLORS["dfcp"]),
        ("mean_emission_iou", COLORS["emission"]),
        ("clade_iou", COLORS["dfcp"]),
        ("emission_clade_iou", COLORS["emission"]),
        ("cluster_purity", COLORS["green"]),
    )
    rows = [(key, color) for key, color in specs if key in result]
    if not rows:
        return None
    fig = go.Figure(go.Bar(
        x=[display_name(key) for key, _ in rows],
        y=[result[key] for key, _ in rows],
        marker_color=[color for _, color in rows],
        text=[f"{result[key]:.3f}" for key, _ in rows],
        textposition="outside",
        hovertemplate="%{x}<br>%{y:.5f}<extra></extra>",
    ))
    fig.update_layout(title="Partition and tree agreement")
    fig.update_yaxes(title_text="similarity / purity", range=[0, 1.08])
    return style_figure(fig)


def timing_figure(result: dict[str, Any]) -> go.Figure | None:
    timings = [
        (key, value) for key, value in result.items()
        if key.startswith("t_") and isinstance(value, (int, float)) and not isinstance(value, bool)
    ]
    train_log = result.get("train_log", [])
    if train_log:
        timings.append(("training total", sum(step.get("t_step", 0) for step in train_log)))
    if not timings:
        return None
    timings.sort(key=lambda item: item[1])
    fig = go.Figure(go.Bar(
        x=[value for _, value in timings],
        y=[display_name(key) for key, _ in timings],
        orientation="h",
        marker_color=COLORS["dfcp"],
        text=[f"{value:,} ms" for _, value in timings],
        textposition="outside",
        hovertemplate="%{y}<br>%{x:,} ms<extra></extra>",
    ))
    fig.update_layout(title="Runtime breakdown")
    fig.update_xaxes(title_text="milliseconds")
    return style_figure(fig, max(340, 80 + 46 * len(timings)))


def kde(values: list[float], grid: list[float]) -> list[float]:
    n = len(values)
    mean = sum(values) / n
    variance = sum((value - mean) ** 2 for value in values) / max(n - 1, 1)
    standard_deviation = math.sqrt(variance)
    bandwidth = max(1.06 * standard_deviation * n ** -0.2, 0.035)
    scale = 1.0 / (n * bandwidth * math.sqrt(2.0 * math.pi))
    return [
        scale * sum(
            math.exp(-0.5 * ((x - value) / bandwidth) ** 2)
            + math.exp(-0.5 * ((x + value) / bandwidth) ** 2)
            for value in values
        )
        for x in grid
    ]


def rgba(color: str, alpha: float) -> str:
    return "rgba({},{},{},{})".format(
        int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16), alpha
    )


def height_figure(result: dict[str, Any]) -> go.Figure | None:
    series = (
        ("DFCP clusters", result.get("dfcp_clade_heights"), COLORS["dfcp"]),
        ("Emissions", result.get("emission_clade_heights"), COLORS["emission"]),
    )
    available = [(name, values, color) for name, values, color in series if values]
    if not available:
        return None
    transformed = [
        (name, [math.log10(1.0 + max(0.0, value)) for value in values], color)
        for name, values, color in available
    ]
    all_values = [value for _, values, _ in transformed for value in values]
    lower = 0.0
    upper = max(all_values)
    padding = max((upper - lower) * 0.025, 0.04)
    upper += padding
    grid = [lower + (upper - lower) * i / 399 for i in range(400)]
    heights = [10.0 ** value - 1.0 for value in grid]
    fig = go.Figure()
    for name, values, color in transformed:
        fig.add_trace(go.Scatter(
            x=[height + 1.0 for height in heights],
            y=kde(values, grid),
            customdata=heights,
            name=name,
            mode="lines",
            fill="tozeroy",
            line={"color": color, "width": 3},
            fillcolor=rgba(color, 0.24),
            hovertemplate=(
                "coalescent height %{customdata:.4g} generations<br>"
                "density %{y:.4f}<extra>" + name + "</extra>"
            ),
        ))
    style_figure(fig, 540)
    fig.update_layout(
        title="Best-matching clade height density",
        hovermode="x",
        legend={"orientation": "h", "x": 1.0, "xanchor": "right", "y": 1.18},
        margin={"l": 58, "r": 24, "t": 92, "b": 50},
    )
    tick_generations = (0, 10, 100, 1_000, 10_000, 100_000)
    tick_values = [1 + value for value in tick_generations if math.log10(1 + value) <= upper]
    tick_text = [
        "0" if value == 0 else f"{value // 1000}k" if value >= 1000 else str(value)
        for value in tick_generations[:len(tick_values)]
    ]
    fig.update_xaxes(
        title_text="coalescent height (generations before present, log scale)",
        type="log", tickmode="array", tickvals=tick_values, ticktext=tick_text,
    )
    fig.update_yaxes(title_text="density")
    return fig


def figure_html(fig: go.Figure | None, div_id: str) -> str:
    if fig is None:
        return ""
    return pio.to_html(
        fig, full_html=False, include_plotlyjs=False,
        config={"responsive": True, "displaylogo": False}, div_id=div_id,
    )


PLOTTED_SCALARS = {
    "viterbi_impute_acc", "fwd_bkwd_impute_acc", "mode_impute_acc",
    "mean_iou", "mean_emission_iou", "clade_iou", "emission_clade_iou",
    "cluster_purity", "t_init", "t_viterbi_impute", "t_fwd_bkwd_impute",
    "t_parsimony", "t_clade_iou", "t_iou",
}


def scalar_table(result: dict[str, Any]) -> str:
    rows = []
    for key, value in result.items():
        if (
            key == "seq_file" or key in PLOTTED_SCALARS
            or not isinstance(value, (int, float)) or isinstance(value, bool)
        ):
            continue
        suffix = " ms" if key.startswith("t_") else ""
        rows.append(
            f"<tr><th>{html.escape(display_name(key))}</th>"
            f"<td>{html.escape(format_value(value) + suffix)}</td></tr>"
        )
    return "".join(rows)


def parameter_cards(result: dict[str, Any]) -> str:
    params = result.get("params", {})
    if not isinstance(params, dict):
        return ""
    cards = []
    for key, value in params.items():
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            cards.append(
                f'<div class="mini-stat"><span>{html.escape(display_name(key))}</span>'
                f"<strong>{html.escape(format_value(value))}</strong></div>"
            )
    return "".join(cards)


SUBGRAPH_RE = re.compile(r"\s*subgraph l(\d+) \{")


def split_dot_loci(path: Path) -> list[tuple[int, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"missing tree DOT file: {path}")
    loci: list[tuple[int, str]] = []
    lines = path.read_text().splitlines()
    line_index = 0
    while line_index < len(lines):
        match = SUBGRAPH_RE.fullmatch(lines[line_index])
        if match is None:
            line_index += 1
            continue
        locus = int(match.group(1))
        body: list[str] = []
        depth = 1
        line_index += 1
        while line_index < len(lines) and depth > 0:
            line = lines[line_index]
            depth += line.count("{") - line.count("}")
            if depth > 0:
                body.append(line)
            line_index += 1
        if depth != 0:
            raise ValueError(f"unterminated locus {locus} subgraph in {path}")
        dot = "\n".join((
            f"digraph locus_{locus} {{",
            '    graph [bgcolor="transparent", pad="0.12", ordering="out"];',
            '    node [style="filled", fontsize="14", penwidth="1.35"];',
            '    edge [color="#4F5F73", penwidth="1.35", arrowsize="0.75"];',
            *body,
            "}",
        ))
        loci.append((locus, dot))
    if not loci:
        raise ValueError(f"{path} contains no DFCP locus subgraphs")
    return loci


def render_dot_svg(dot: str, locus: int) -> str:
    try:
        result = subprocess.run(
            ["dot", "-Tsvg"], input=dot, capture_output=True, text=True, check=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError("Graphviz 'dot' is required for --tree-dot") from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"Graphviz failed for locus {locus}: {error.stderr.strip()}") from error
    svg = result.stdout[result.stdout.find("<svg"):]
    if not svg:
        raise RuntimeError(f"Graphviz emitted no SVG for locus {locus}")
    prefix = f"tree-locus-{locus}-"
    svg = re.sub(r'\bid="([^"]+)"', lambda match: f'id="{prefix}{match.group(1)}"', svg)
    svg = re.sub(
        r'((?:xlink:)?href)="#([^"]+)"',
        lambda match: f'{match.group(1)}="#{prefix}{match.group(2)}"',
        svg,
    )
    return svg.replace("<svg ", f'<svg class="tree-svg" data-locus="{locus}" ', 1)


def tree_report_html(path: Path | None) -> str:
    if path is None:
        return ""
    loci = split_dot_loci(path)
    frames = []
    for index, (locus, dot) in enumerate(loci):
        hidden = "" if index == 0 else " hidden"
        frames.append(
            f'<div class="tree-frame" data-index="{index}" data-locus="{locus}"{hidden}>'
            f"{render_dot_svg(dot, locus)}</div>"
        )
    locus_labels = json.dumps([locus for locus, _ in loci])
    return f"""
<section class="plot-section wide tree-section">
  <div class="tree-heading"><h2>Marginal trees</h2>
    <div class="tree-zoom"><button type="button" data-tree-action="minus">−</button>
      <button type="button" data-tree-action="readable">Readable</button>
      <button type="button" data-tree-action="fit">Fit</button>
      <button type="button" data-tree-action="plus">+</button></div>
  </div>
  <div class="tree-canvas">{''.join(frames)}</div>
  <div class="tree-slider-row"><input id="tree-locus-slider" type="range" min="0" max="{len(loci)-1}" value="0" step="1">
    <output id="tree-locus-value">locus {loci[0][0]}</output></div>
</section>
<script>
(() => {{
  const loci = {locus_labels};
  const frames = [...document.querySelectorAll('.tree-frame')];
  const canvas = document.querySelector('.tree-canvas');
  const slider = document.getElementById('tree-locus-slider');
  const output = document.getElementById('tree-locus-value');
  let treeWidth = 3600;
  function activeSvg() {{ return frames[Number(slider.value)].querySelector('svg'); }}
  function applyWidth() {{ activeSvg().style.width = `${{treeWidth}}px`; }}
  function centerTree() {{ canvas.scrollLeft = Math.max(0, (canvas.scrollWidth - canvas.clientWidth) / 2); }}
  function showTree() {{
    frames.forEach((frame, index) => frame.hidden = index !== Number(slider.value));
    output.value = `locus ${{loci[Number(slider.value)]}}`;
    applyWidth();
    requestAnimationFrame(centerTree);
  }}
  slider.addEventListener('input', showTree);
  document.querySelectorAll('[data-tree-action]').forEach(button => button.addEventListener('click', () => {{
    const action = button.dataset.treeAction;
    if (action === 'minus') treeWidth = Math.max(800, treeWidth / 1.25);
    if (action === 'plus') treeWidth = Math.min(8000, treeWidth * 1.25);
    if (action === 'readable') treeWidth = 3600;
    if (action === 'fit') treeWidth = Math.max(700, canvas.clientWidth - 24);
    applyWidth();
    requestAnimationFrame(centerTree);
  }}));
  showTree();
}})();
</script>"""


CSS = r"""
:root { color-scheme:light; --ink:#2a3f5f; --muted:#697386; --line:#e5e7eb; }
* { box-sizing:border-box; }
body { margin:0; color:var(--ink); background:white;
  font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
.run-sidebar { position:fixed; inset:0 auto 0 0; z-index:5; width:310px; overflow-y:auto;
  padding:24px 20px; border-right:1px solid var(--line); background:#fafbfc; }
.run-sidebar h1 { margin:0 0 8px; font-size:19px; font-weight:650; }
.seq-file { margin:0; color:var(--muted); font:12px/1.5 ui-monospace,monospace; overflow-wrap:anywhere; }
.sidebar-section { margin-top:24px; }
.sidebar-section h2 { margin:0 0 9px; font-size:14px; font-weight:650; }
.sidebar-note { margin:0; color:var(--muted); font-size:12px; line-height:1.5; }
.sidebar-table { width:100%; border-collapse:collapse; font-size:12px; }
.sidebar-table th,.sidebar-table td { padding:6px 2px; border-bottom:1px solid var(--line); }
.sidebar-table th { text-align:left; color:var(--muted); font-weight:500; }
.sidebar-table td { text-align:right; font:600 12px ui-monospace,monospace; }
.sidebar-params .mini-stat { display:grid; grid-template-columns:minmax(0,1fr) auto; gap:8px;
  padding:7px 2px; border-bottom:1px solid var(--line); font-size:12px; }
.sidebar-params .mini-stat span { color:var(--muted); }
.sidebar-params .mini-stat strong { font:600 12px ui-monospace,monospace; }
.report-main { width:calc(100% - 310px); margin-left:310px; padding:12px 28px 60px; }
.grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); column-gap:28px; row-gap:34px; }
.plot-section { min-width:0; }
.plot-section.wide { grid-column:1/-1; }
.tree-section { margin-top:44px; }
.tree-heading { display:flex; align-items:flex-end; justify-content:space-between; gap:18px; margin:0 8px 12px; }
.tree-heading h2 { margin:0; font-size:17px; }
.tree-zoom { display:flex; gap:5px; flex-shrink:0; }
.tree-zoom button { padding:6px 10px; border:1px solid #aeb5c0; border-radius:4px;
  color:var(--ink); background:#f4f5f7; cursor:pointer; }
.tree-zoom button:hover { border-color:#3366cc; color:#3366cc; }
.tree-canvas { height:720px; overflow:auto; border:1px solid var(--line); background:white; }
.tree-frame { min-width:100%; padding:10px 12px 16px; }
.tree-frame[hidden] { display:none; }
.tree-svg { display:block; height:auto; max-width:none; margin:0 auto; }
.tree-slider-row { display:grid; grid-template-columns:minmax(0,1fr) 72px; align-items:center;
  gap:14px; padding:14px 8px 0; }
.tree-slider-row input { width:100%; }
.tree-slider-row output { font:600 11px ui-monospace,monospace; }
@media (max-width:950px) {
  .run-sidebar { width:260px; }
  .report-main { width:calc(100% - 260px); margin-left:260px; padding-inline:18px; }
  .grid { grid-template-columns:1fr; }
  .plot-section.wide { grid-column:auto; }
}
"""


def write_report(result: dict[str, Any], output: Path, tree_dot: Path | None) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    plotly_asset = ensure_plotly_asset(output)
    plots = {
        "elbo": figure_html(elbo_figure(result), "dfcp-elbo"),
        "iteration-time": figure_html(iteration_time_figure(result), "dfcp-iteration-time"),
        "parameters": figure_html(parameter_figure(result), "dfcp-parameters"),
        "imputation": figure_html(imputation_figure(result), "dfcp-imputation"),
        "metrics": figure_html(metric_figure(result), "dfcp-metrics"),
        "timing": figure_html(timing_figure(result), "dfcp-timing"),
        "heights": figure_html(height_figure(result), "dfcp-heights"),
    }
    panels = []
    for key in ("elbo", "iteration-time"):
        if plots[key]:
            panels.append(f'<section class="plot-section">{plots[key]}</section>')
    if plots["parameters"]:
        panels.append(f'<section class="plot-section">{plots["parameters"]}</section>')
    for key in ("imputation", "metrics", "timing"):
        if plots[key]:
            panels.append(f'<section class="plot-section">{plots[key]}</section>')
    if plots["heights"]:
        panels.append(f'<section class="plot-section wide">{plots["heights"]}</section>')

    scalars = scalar_table(result)
    scalar_sidebar = (
        '<section class="sidebar-section"><h2>Additional metrics</h2>'
        f'<table class="sidebar-table">{scalars}</table></section>' if scalars else ""
    )
    params = parameter_cards(result)
    parameter_sidebar = (
        '<section class="sidebar-section"><h2>Scalar parameters</h2>'
        f'<div class="sidebar-params">{params}</div></section>' if params else ""
    )
    tree_sidebar = (
        '<section class="sidebar-section"><h2>Tree display</h2>'
        '<p class="sidebar-note">Node labels are tree node IDs. Leaf color is the DFCP cluster; '
        'leaf shape is the emission. Use Readable for labeled inspection or Fit for an overview.</p>'
        '</section>' if tree_dot is not None else ""
    )
    tree_report = tree_report_html(tree_dot)

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DFCP run report</title><script src="{html.escape(plotly_asset, quote=True)}"></script><style>{CSS}</style></head>
<body><aside class="run-sidebar"><h1>DFCP run</h1>
<p class="seq-file">{html.escape(result['seq_file'])}</p>{parameter_sidebar}{scalar_sidebar}{tree_sidebar}</aside>
<main class="report-main"><div class="grid">{''.join(panels)}</div>{tree_report}</main></body></html>
"""
    output.write_text(page)


def main() -> None:
    args = parse_args()
    result = load_result(args.json)
    write_report(result, args.output, args.tree_dot)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
