import json
import re
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dfcp import build_dfcp, run_dfcp, get_dfcp_parser
from plotly_html import ensure_plotly_asset


ERROR_RE = re.compile(r"\.txt\.gz_([0-9.]+)_([0-9.]+)\.txt\.gz_")


def get_seq_label(seq_file: Path) -> str:
    match = ERROR_RE.search(seq_file.name)
    if match is None:
        return "Baseline (no injected errors)"
    return f"Bit flip {match.group(1)}, switch {match.group(2)}"


def get_seq_sort_key(seq_file: Path) -> tuple[bool, float, float]:
    match = ERROR_RE.search(seq_file.name)
    if match is None:
        return False, 0.0, 0.0
    return True, float(match.group(1)), float(match.group(2))


if __name__ == "__main__":
    p = get_dfcp_parser()
    p.set_defaults(seq_file=None)
    p.add_argument(
        "--seq_dir", type=Path,
        default=Path("data/examples/simulated/SIM1_LEN500_NHAPS100"),
    )
    p.add_argument("--seq_files", type=Path, nargs="+")
    p.add_argument(
        "--load", type=Path, metavar="JSON",
        help="skip DFCP runs and rebuild the Plotly report from saved results",
    )
    p.add_argument("--json", type=Path, default=Path("output/init.json"))
    p.add_argument("--output", type=Path, default=Path("docs/init.html"))
    args = p.parse_args()

    INIT_METHODS = (
        "block",
        ("pbwt", 2), ("pbwt", 5), ("pbwt", 10), ("pbwt", 20), ("pbwt", 50), ("pbwt", 100),
        "viterbi"
    )
    PBWT_MATCH_LENGTHS = tuple(i[1] for i in INIT_METHODS if i[0] == "pbwt")
    PBWT_COLORS = dict(zip(
        PBWT_MATCH_LENGTHS,
        ("#F6C945", "#F39C34", "#E7682B", "#C9342C", "#8E1B1B", "#5E0B0B"),
    ))
    LINE_STYLES = {
        "block": {"color": "#3366CC", "width": 3, "dash": "dash"},
        "viterbi": {"color": "#2A9D8F", "width": 3, "dash": "dot"},
    }
    SHARED_METRICS = (
        "t_init",
        "dfcp_impute_acc",
        "clade_iou", "mean_excess_parsimony",
        "mean_iou",
        "mean_clusters",
    )
    POST_METRICS = ("elbo", "iterations")
    METRICS = (
        "dfcp_impute_acc", "elbo",
        "clade_iou", "mean_excess_parsimony",
        "mean_iou", "mean_clusters",
        "t_init", "iterations",
    )
    PLOT_TITLES = {
        "elbo": "final elbo (post-training)",
        "iterations": "# iterations to convergence (post-training)",
    }
    PHASES = {"pre": 1, "post": 0}
    MODES = {
        "hard": {"noisy": 0, "soft": 0},
        "noisy": {"noisy": 1, "soft": 0},
        "soft": {"noisy": 0, "soft": 1},
    }
    MASKS = (0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)

    def empty_data() -> dict:
        return {
            seq_idx: {
                mode: {
                    phase: {i: {m: [] for m in METRICS} for i in INIT_METHODS}
                    for phase in PHASES
                }
                for mode in MODES
            }
            for seq_idx in range(len(seq_files))
        }

    def method_id(init_method: str | tuple[str, int]) -> str:
        return (
            f"pbwt:{init_method[1]}"
            if isinstance(init_method, tuple) else init_method
        )

    def method_from_id(identifier: str) -> str | tuple[str, int]:
        if identifier.startswith("pbwt:"):
            return "pbwt", int(identifier.split(":", 1)[1])
        if identifier not in ("block", "viterbi"):
            raise ValueError(f"unknown initialization method in JSON: {identifier}")
        return identifier

    if args.load is not None:
        payload = json.loads(args.load.read_text())
        if payload.get("version") != 1 or not isinstance(payload.get("series"), list):
            raise ValueError(f"{args.load} is not a DFCP initialization result")
        seq_files = [Path(path) for path in payload["config"]["seq_files"]]
        if not seq_files:
            raise ValueError(f"{args.load} contains no sequence files")
        MASKS = tuple(float(mask) for mask in payload["config"]["masks"])
        if not MASKS:
            raise ValueError(f"{args.load} contains no mask values")
        data = empty_data()
        for series in payload["series"]:
            init_method = method_from_id(series["method"])
            data[series["dataset"]][series["mode"]][series["phase"]][init_method] = (
                series["metrics"]
            )
        print(f"loaded {args.load}")
    else:
        if args.seq_files is not None:
            seq_files = args.seq_files
        elif args.seq_file is not None:
            seq_files = [Path(args.seq_file)]
        else:
            seq_files = sorted(args.seq_dir.glob("haps_*.txt"), key=get_seq_sort_key)
        if not seq_files:
            raise RuntimeError("no haplotype sequence files found")

        data = empty_data()
        build_dfcp()
        for seq_idx, seq_file in enumerate(seq_files):
            for mode, mode_flags in MODES.items():
                for init_method in INIT_METHODS:
                    print(f"{seq_file} {mode} {init_method}")
                    for mask in MASKS:
                        for phase, init_only in PHASES.items():
                            res = run_dfcp(
                                str(seq_file),

                                val=0.2, mask=mask,

                                tree=args.tree, variant_pos_fname=args.variant_pos_fname,
                                variant_start_pos=args.variant_start_pos,
                                clade_beta=args.clade_beta,

                                noisy=mode_flags["noisy"], soft=mode_flags["soft"],
                                block_init=int(init_method == "block"),
                                pbwt_init=int(init_method[0] == "pbwt"),
                                pbwt_match_len=(
                                    init_method[1] if init_method[0] == "pbwt" else None
                                ),
                                init_only=init_only,
                            )
                            for metric in SHARED_METRICS:
                                data[seq_idx][mode][phase][init_method][metric].append(
                                    res[metric]
                                )
                            if phase == "post":
                                train_log = res["train_log"]
                                data[seq_idx][mode][phase][init_method]["elbo"].append(
                                    train_log[-1]["elbo"]
                                )
                                data[seq_idx][mode][phase][init_method]["iterations"].append(
                                    len(train_log)
                                )

        payload = {
            "version": 1,
            "config": {
                "seq_files": [str(path) for path in seq_files],
                "masks": list(MASKS),
                "val": 0.2,
                "tree": args.tree,
                "variant_pos_fname": args.variant_pos_fname,
                "variant_start_pos": args.variant_start_pos,
                "clade_beta": args.clade_beta,
            },
            "series": [
                {
                    "dataset": seq_idx,
                    "mode": mode,
                    "phase": phase,
                    "method": method_id(init_method),
                    "metrics": data[seq_idx][mode][phase][init_method],
                }
                for seq_idx in range(len(seq_files))
                for mode in MODES
                for phase in PHASES
                for init_method in INIT_METHODS
            ],
        }
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
        print(f"wrote {args.json}")

    ncols = 2
    nrows = (len(METRICS) + ncols - 1) // ncols
    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=tuple(PLOT_TITLES.get(metric, metric) for metric in METRICS),
        vertical_spacing=0.10, horizontal_spacing=0.08,
    )
    for plot_idx, metric in enumerate(METRICS):
        row = plot_idx // ncols + 1
        col = plot_idx % ncols + 1
        phases = ("post",) if metric in POST_METRICS else PHASES
        for seq_idx, seq_file in enumerate(seq_files):
            for mode in MODES:
                for phase in phases:
                    trace_phase = "post_only" if metric in POST_METRICS else phase
                    for init_method in INIT_METHODS:
                        is_pbwt = init_method[0] == "pbwt"
                        method = "pbwt" if is_pbwt else init_method
                        name = f"PBWT {init_method[1]}" if is_pbwt else init_method
                        line = (
                            {"color": PBWT_COLORS[init_method[1]], "width": 2}
                            if is_pbwt else LINE_STYLES[init_method]
                        )
                        fig.add_trace(
                            go.Scatter(
                                x=MASKS,
                                y=data[seq_idx][mode][phase][init_method][metric],
                                name=name,
                                legendgroup=name,
                                showlegend=False,
                                visible=(
                                    seq_idx == 0 and mode == "hard" and
                                    trace_phase in ("pre", "post_only")
                                ),
                                line=line,
                                meta={
                                    "dataset": seq_idx,
                                    "dataset_file": seq_file.name,
                                    "dataset_label": get_seq_label(seq_file),
                                    "method": method,
                                    "match_len": init_method[1] if is_pbwt else None,
                                    "mode": mode,
                                    "phase": trace_phase,
                                },
                            ),
                            row=row,
                            col=col,
                        )

    fig.update_xaxes(title_text="mask frac")
    fig.update_layout(
        height=nrows * 450,
        title=(
            "DFCP init methods: pre/post-training metrics"
        ),
        margin={"b": 60, "t": 100},
    )

    post_script = r"""
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
const pbwtLegend = pbwtMatchLengths.map(matchLength => {
    const trace = plot.data.find(candidate =>
        candidate.meta.method === 'pbwt' && candidate.meta.match_len === matchLength
    );
    return legendEntry(`PBWT ${matchLength}`, trace.line.color);
}).join('');
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
window.addEventListener('resize', () => Plotly.Plots.resize(plot));
"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        args.output,
        include_plotlyjs=ensure_plotly_asset(args.output),
        config={"responsive": True},
        post_script=post_script,
    )
    print(f"wrote {args.output}")
