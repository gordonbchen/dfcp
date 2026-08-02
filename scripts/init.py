import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dfcp import build_dfcp, run_dfcp, get_dfcp_parser


if __name__ == "__main__":
    p = get_dfcp_parser()
    args = p.parse_args()

    build_dfcp()

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
        "t_init", "iterations",
        "dfcp_impute_acc", "elbo",
        "clade_iou", "mean_excess_parsimony",
        "mean_iou", "mean_clusters",
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
    data = {
        mode: {
            phase: {i: {m: [] for m in METRICS} for i in INIT_METHODS}
            for phase in PHASES
        }
        for mode in MODES
    }

    MASKS = (0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99)

    for mode, mode_flags in MODES.items():
        for init_method in INIT_METHODS:
            for mask in MASKS:
                for phase, init_only in PHASES.items():
                    res = run_dfcp(
                        args.seq_file,

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
                        data[mode][phase][init_method][metric].append(res[metric])
                    if phase == "post":
                        train_log = res["train_log"]
                        data[mode][phase][init_method]["elbo"].append(
                            train_log[-1]["elbo"]
                        )
                        data[mode][phase][init_method]["iterations"].append(
                            len(train_log)
                        )

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
                            y=data[mode][phase][init_method][metric],
                            name=name,
                            legendgroup=name,
                            showlegend=False,
                            visible=(
                                mode == "hard" and
                                trace_phase in ("pre", "post_only")
                            ),
                            line=line,
                            meta={
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
        title=f"DFCP init methods: pre/post-training metrics<br>seq_file: {args.seq_file}",
        margin={"b": 60, "t": 100},
    )

    post_script = r"""
const plot = document.getElementById('{plot_id}');
const pbwtMatchLengths = [...new Set(
    plot.data
        .filter(trace => trace.meta.method === 'pbwt')
        .map(trace => trace.meta.match_len)
)].sort((a, b) => a - b);

let selectedMatchLength = pbwtMatchLengths[0];
let showAllPbwt = true;
let trainingPhase = 'pre';
let modelMode = 'hard';

function updateVisibleTraces() {
    const visible = plot.data.map(trace =>
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

const sidebar = document.createElement('aside');
sidebar.id = 'dfcp-sidebar';
sidebar.innerHTML = `
    <h2>Controls</h2>
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
    fig.write_html("init.html", config={"responsive": True}, post_script=post_script)
