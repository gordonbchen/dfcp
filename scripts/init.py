"""Run and visualize the DFCP initialization-method comparison."""

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dfcp import build_dfcp, run_dfcp
from plotly_html import ensure_plotly_asset
from seq_file_name import get_seq_label, get_seq_sort_key


DEFAULT_SEQ_DIR = Path("data/examples/simulated/SIM1_LEN500_NHAPS100")
MASKS = (0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
MODES = ("hard", "noisy", "soft")
PHASES = {"pre": 1, "post": 0}
SHARED_METRICS = (
    "t_init",
    "viterbi_impute_acc",
    "clade_iou",
    "mean_excess_parsimony",
    "mean_iou",
    "mean_clusters",
)
POST_METRICS = ("elbo", "iterations")
METRICS = (
    "viterbi_impute_acc",
    "elbo",
    "clade_iou",
    "mean_excess_parsimony",
    "mean_iou",
    "mean_clusters",
    "t_init",
    "iterations",
)
PLOT_TITLES = {
    "elbo": "final elbo (post-training)",
    "iterations": "# iterations to convergence (post-training)",
}
PBWT_COLORS = {
    5: "#F39C34",
    10: "#E7682B",
    20: "#C9342C",
    50: "#8E1B1B",
    100: "#5E0B0B",
}
LINE_STYLES = {
    "block": {"color": "#3366CC", "width": 3, "dash": "dash"},
    "viterbi": {"color": "#2A9D8F", "width": 3, "dash": "dot"},
}


@dataclass(frozen=True)
class InitMethod:
    name: str
    match_len: int | None = None
    match_curr: bool | None = None

    @property
    def label(self) -> str:
        if self.name == "pbwt":
            suffix = " (allow curr mismatch)" if not self.match_curr else ""
            return f"PBWT {self.match_len}{suffix}"
        return self.name


INIT_METHODS = (
    InitMethod("block"),
    InitMethod("pbwt", 5, True),
    InitMethod("pbwt", 10, True),
    InitMethod("pbwt", 20, True),
    InitMethod("pbwt", 20, False),
    InitMethod("pbwt", 50, True),
    InitMethod("pbwt", 100, True),
    InitMethod("viterbi"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare DFCP initialization methods and save reusable results."
    )
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--seq_file", type=Path)
    inputs.add_argument("--seq_files", type=Path, nargs="+")
    inputs.add_argument("--seq_dir", type=Path)
    inputs.add_argument(
        "--load", type=Path, metavar="JSON",
        help="skip DFCP and rebuild the Plotly report from saved results",
    )
    parser.add_argument("--val", type=float, default=0.2)
    parser.add_argument("--tree", type=Path)
    parser.add_argument("--variant_pos_fname", type=Path)
    parser.add_argument("--variant_start_pos", type=int)
    parser.add_argument("--clade_beta", type=float)
    parser.add_argument("--json", type=Path, default=Path("output/init.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/init.html"))
    return parser.parse_args()


def resolve_seq_files(args: argparse.Namespace) -> list[Path]:
    if args.seq_files:
        files = args.seq_files
    elif args.seq_file:
        files = [args.seq_file]
    else:
        seq_dir = args.seq_dir or DEFAULT_SEQ_DIR
        files = sorted(seq_dir.glob("haps_*.txt*"), key=get_seq_sort_key)
    if not files:
        raise ValueError("no haplotype sequence files were found")
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(f"missing sequence file: {path}")
    return files


def collect_experiment(args: argparse.Namespace, seq_files: list[Path]) -> dict:
    datasets = []
    for seq_file in seq_files:
        dataset = {
            "seq_file": str(seq_file),
            "series": [],
        }
        for mode in MODES:
            for method in INIT_METHODS:
                print(f"{seq_file} {mode} {method.label}", flush=True)
                phase_series = {
                    phase: {
                        "mode": mode,
                        "phase": phase,
                        "method": method.name,
                        "match_len": method.match_len,
                        "match_curr": method.match_curr,
                        "metrics": {metric: [] for metric in METRICS},
                    }
                    for phase in PHASES
                }
                for mask in MASKS:
                    for phase, init_only in PHASES.items():
                        result = run_dfcp(
                            seq_file,
                            retries=1,
                            val=args.val,
                            mask=mask,
                            tree=args.tree,
                            variant_pos_fname=args.variant_pos_fname,
                            variant_start_pos=args.variant_start_pos,
                            clade_beta=args.clade_beta,
                            mode=mode,
                            init=method.name,
                            pbwt_match_len=method.match_len,
                            pbwt_match_curr=(
                                int(method.match_curr) if method.name == "pbwt" else None
                            ),
                            init_only=init_only,
                        )
                        metrics = phase_series[phase]["metrics"]
                        for metric in SHARED_METRICS:
                            metrics[metric].append(result[metric])
                        if phase == "post":
                            train_log = result["train_log"]
                            metrics["elbo"].append(train_log[-1]["elbo"])
                            metrics["iterations"].append(len(train_log))
                dataset["series"].extend(phase_series.values())
        datasets.append(dataset)

    return {
        "config": {
            "masks": list(MASKS),
            "val": args.val,
            "tree": str(args.tree) if args.tree else None,
            "variant_pos_fname": (
                str(args.variant_pos_fname) if args.variant_pos_fname else None
            ),
            "variant_start_pos": args.variant_start_pos,
            "clade_beta": args.clade_beta,
        },
        "datasets": datasets,
    }


def load_experiment(path: Path) -> dict:
    experiment = json.loads(path.read_text())
    datasets = experiment.get("datasets")
    masks = experiment.get("config", {}).get("masks")
    if not isinstance(datasets, list) or not datasets or not isinstance(masks, list):
        raise ValueError(f"{path} is not a DFCP initialization result")
    return experiment


def write_json(experiment: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(experiment, indent=0, allow_nan=False) + "\n")


def make_figure(experiment: dict) -> go.Figure:
    masks = experiment["config"]["masks"]
    datasets = sorted(
        experiment["datasets"],
        key=lambda dataset: get_seq_sort_key(Path(dataset["seq_file"])),
    )
    ncols = 2
    nrows = math.ceil(len(METRICS) / ncols)
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        subplot_titles=tuple(PLOT_TITLES.get(metric, metric) for metric in METRICS),
        vertical_spacing=0.10,
        horizontal_spacing=0.08,
    )
    for plot_index, metric in enumerate(METRICS):
        row = plot_index // ncols + 1
        col = plot_index % ncols + 1
        for dataset_index, dataset in enumerate(datasets):
            seq_file = Path(dataset["seq_file"])
            for series in dataset["series"]:
                if metric in POST_METRICS and series["phase"] != "post":
                    continue
                trace_phase = (
                    "post_only" if metric in POST_METRICS else series["phase"]
                )
                is_pbwt = series["method"] == "pbwt"
                match_curr = series.get("match_curr", True)
                line = (
                    {
                        "color": PBWT_COLORS[series["match_len"]],
                        "width": 2,
                        "dash": "solid" if match_curr else "dash",
                    }
                    if is_pbwt else LINE_STYLES[series["method"]]
                )
                fig.add_trace(
                    go.Scatter(
                        x=masks,
                        y=series["metrics"][metric],
                        name=(
                            f"PBWT {series['match_len']}"
                            + ("" if match_curr else " (allow curr mismatch)")
                            if is_pbwt else series["method"]
                        ),
                        showlegend=False,
                        visible=(
                            dataset_index == 0
                            and series["mode"] == "hard"
                            and trace_phase in ("pre", "post_only")
                        ),
                        line=line,
                        meta={
                            "dataset": dataset_index,
                            "dataset_file": seq_file.name,
                            "dataset_label": get_seq_label(seq_file),
                            "method": series["method"],
                            "match_len": series["match_len"],
                            "match_curr": match_curr if is_pbwt else None,
                            "mode": series["mode"],
                            "phase": trace_phase,
                        },
                    ),
                    row=row,
                    col=col,
                )
    fig.update_xaxes(title_text="mask frac")
    fig.update_layout(
        height=nrows * 450,
        title="DFCP init methods: pre/post-training metrics",
        margin={"b": 60, "t": 100},
    )
    return fig


POST_SCRIPT = r"""
const plot = document.getElementById('{plot_id}');
const datasets = [...new Map(plot.data.map(trace => [
    trace.meta.dataset,
    {
        id: trace.meta.dataset,
        label: trace.meta.dataset_label,
        file: trace.meta.dataset_file,
    },
])).values()];
const pbwtMatchLengths = [...new Set(
    plot.data
        .filter(trace => trace.meta.method === 'pbwt')
        .map(trace => trace.meta.match_len)
)].sort((a, b) => a - b);

let selectedMatchLength = pbwtMatchLengths[0];
let showAllPbwt = true;
let trainingPhase = 'pre';
let modelMode = 'hard';
let selectedDataset = datasets[0].id;

function updateVisibleTraces() {
    const visible = plot.data.map(trace =>
        trace.meta.dataset === selectedDataset &&
        trace.meta.mode === modelMode &&
        (trace.meta.phase === trainingPhase || trace.meta.phase === 'post_only') &&
        (trace.meta.method !== 'pbwt' || showAllPbwt ||
            trace.meta.match_len === selectedMatchLength)
    );
    Plotly.restyle(plot, {visible: visible});
}

function bindButtons(containerId, onChange) {
    const buttons = document.querySelectorAll(`#${containerId} button`);
    buttons.forEach(button => {
        button.addEventListener('click', () => {
            buttons.forEach(other => other.classList.remove('active'));
            button.classList.add('active');
            onChange(button.dataset.value);
            updateVisibleTraces();
        });
    });
}

function legendEntry(label, color, dash = 'solid', width = 2) {
    const borderStyle = dash === 'dash' ? 'dashed' : dash === 'dot' ? 'dotted' : 'solid';
    return `<div class="dfcp-legend-entry">
        <span style="border-top:${width}px ${borderStyle} ${color}"></span>
        <span>${label}</span>
    </div>`;
}

const blockTrace = plot.data.find(trace => trace.meta.method === 'block');
const viterbiTrace = plot.data.find(trace => trace.meta.method === 'viterbi');
const pbwtTraces = [...new Map(
    plot.data
        .filter(trace => trace.meta.method === 'pbwt')
        .map(trace => [`${trace.meta.match_len}:${trace.meta.match_curr}`, trace])
).values()];
const pbwtLegend = pbwtTraces.map(trace => legendEntry(
    trace.name, trace.line.color, trace.line.dash, trace.line.width
)).join('');
const legend = [
    legendEntry('block', blockTrace.line.color, blockTrace.line.dash, blockTrace.line.width),
    pbwtLegend,
    legendEntry(
        'viterbi', viterbiTrace.line.color, viterbiTrace.line.dash, viterbiTrace.line.width
    ),
].join('');
const sliderTicks = pbwtMatchLengths.map(length => `<span>${length}</span>`).join('');
const datasetOptions = datasets.map(dataset =>
    `<option value="${dataset.id}">${dataset.label}</option>`
).join('');

const sidebar = document.createElement('aside');
sidebar.id = 'dfcp-sidebar';
sidebar.innerHTML = `
    <h2>Controls</h2>
    <section>
        <h3>Haplotype dataset</h3>
        <select id="dataset-select">${datasetOptions}</select>
        <div class="dfcp-dataset-file" id="dataset-file">${datasets[0].file}</div>
    </section>
    <section>
        <h3>Model</h3>
        <div class="dfcp-buttons" id="model-controls">
            <button class="active" data-value="hard">Hard</button>
            <button data-value="noisy">Noisy</button>
            <button data-value="soft">Soft</button>
        </div>
    </section>
    <section>
        <h3>Training stage</h3>
        <div class="dfcp-buttons" id="phase-controls">
            <button class="active" data-value="pre">Pre</button>
            <button data-value="post">Post</button>
        </div>
    </section>
    <section>
        <h3>PBWT traces</h3>
        <div class="dfcp-buttons" id="pbwt-display-controls">
            <button class="active" data-value="all">All</button>
            <button data-value="selected">Selected</button>
        </div>
    </section>
    <section>
        <h3>PBWT match length: <output id="pbwt-match-value">${selectedMatchLength}</output></h3>
        <input id="pbwt-match-slider" type="range" min="0"
            max="${pbwtMatchLengths.length - 1}" value="0" step="1">
        <div class="dfcp-slider-ticks">${sliderTicks}</div>
    </section>
    <section>
        <h3>Legend</h3>
        <div class="dfcp-legend">${legend}</div>
    </section>
`;
document.body.prepend(sidebar);

const style = document.createElement('style');
style.textContent = `
    html, body { margin: 0; }
    #dfcp-sidebar {
        position: fixed;
        inset: 0 auto 0 0;
        z-index: 1000;
        box-sizing: border-box;
        width: 280px;
        padding: 18px;
        overflow-y: auto;
        border-right: 1px solid #d9dce1;
        background: rgba(255, 255, 255, 0.98);
        color: #2a3f5f;
        font-family: Arial, sans-serif;
        box-shadow: 2px 0 8px rgba(0, 0, 0, 0.08);
    }
    #dfcp-sidebar h2 { margin: 0 0 18px; font-size: 21px; }
    #dfcp-sidebar h3 { margin: 0 0 8px; font-size: 14px; }
    #dfcp-sidebar section { margin-bottom: 20px; }
    #dataset-select {
        box-sizing: border-box;
        width: 100%;
        padding: 7px 5px;
        border: 1px solid #aeb5c0;
        border-radius: 4px;
        background: white;
        color: #2a3f5f;
    }
    .dfcp-dataset-file {
        margin-top: 6px;
        overflow-wrap: anywhere;
        color: #697386;
        font-size: 11px;
        line-height: 1.3;
    }
    .dfcp-buttons { display: flex; gap: 5px; }
    .dfcp-buttons button {
        flex: 1;
        padding: 7px 5px;
        border: 1px solid #aeb5c0;
        border-radius: 4px;
        background: #f4f5f7;
        color: #2a3f5f;
        cursor: pointer;
    }
    .dfcp-buttons button:hover { background: #e5e9ef; }
    .dfcp-buttons button.active {
        border-color: #3366cc;
        background: #3366cc;
        color: white;
    }
    #pbwt-match-slider { width: 100%; accent-color: #c9342c; }
    .dfcp-slider-ticks {
        display: flex;
        justify-content: space-between;
        margin: 2px 2px 0;
        color: #697386;
        font-size: 11px;
    }
    .dfcp-legend { display: grid; gap: 8px; }
    .dfcp-legend-entry {
        display: grid;
        grid-template-columns: 40px 1fr;
        align-items: center;
        gap: 8px;
        font-size: 13px;
    }
    .dfcp-legend-entry > span:first-child { width: 36px; }
    .dfcp-plot {
        width: calc(100% - 280px) !important;
        margin-left: 280px;
    }
    @media (max-width: 800px) {
        #dfcp-sidebar { width: 230px; padding: 12px; }
        .dfcp-plot { width: calc(100% - 230px) !important; margin-left: 230px; }
    }
`;
document.head.appendChild(style);
plot.classList.add('dfcp-plot');

bindButtons('model-controls', value => { modelMode = value; });
bindButtons('phase-controls', value => { trainingPhase = value; });
bindButtons('pbwt-display-controls', value => { showAllPbwt = value === 'all'; });

const datasetSelect = document.getElementById('dataset-select');
datasetSelect.addEventListener('change', () => {
    selectedDataset = Number(datasetSelect.value);
    document.getElementById('dataset-file').textContent = datasets.find(
        dataset => dataset.id === selectedDataset
    ).file;
    updateVisibleTraces();
});

const slider = document.getElementById('pbwt-match-slider');
slider.addEventListener('input', () => {
    selectedMatchLength = pbwtMatchLengths[Number(slider.value)];
    document.getElementById('pbwt-match-value').value = selectedMatchLength;
    showAllPbwt = false;
    document.querySelectorAll('#pbwt-display-controls button').forEach(button => {
        button.classList.toggle('active', button.dataset.value === 'selected');
    });
    updateVisibleTraces();
});

requestAnimationFrame(() => Plotly.Plots.resize(plot));
window.addEventListener('resize', () => Plotly.Plots.resize(plot));"""


def write_report(experiment: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    make_figure(experiment).write_html(
        output,
        include_plotlyjs=ensure_plotly_asset(output),
        config={"responsive": True},
        post_script=POST_SCRIPT,
    )


def main() -> None:
    args = parse_args()
    if args.load:
        experiment = load_experiment(args.load)
        print(f"loaded {args.load}")
    else:
        if not 0.0 < args.val < 1.0:
            raise ValueError("val must be in (0,1)")
        seq_files = resolve_seq_files(args)
        build_dfcp()
        experiment = collect_experiment(args, seq_files)
        write_json(experiment, args.json)
        print(f"wrote {args.json}")

    write_report(experiment, args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
