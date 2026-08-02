import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dfcp import build_dfcp, run_dfcp, get_dfcp_parser


if __name__ == "__main__":
    p = get_dfcp_parser()
    args = p.parse_args()

    build_dfcp()

    INIT_METHODS = (
        "block",
        ("pbwt", 2), ("pbwt", 5), ("pbwt", 10), ("pbwt", 20), ("pbwt", 50),
        "viterbi"
    )
    PBWT_MATCH_LENGTHS = tuple(i[1] for i in INIT_METHODS if i[0] == "pbwt")
    PBWT_COLORS = dict(zip(
        PBWT_MATCH_LENGTHS,
        ("#F6C945", "#F39C34", "#E7682B", "#C9342C", "#8E1B1B"),
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
    data = {
        phase: {i: {m: [] for m in METRICS} for i in INIT_METHODS}
        for phase in PHASES
    }

    MASKS = (0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99)

    for init_method in INIT_METHODS:
        for mask in MASKS:
            for phase, init_only in PHASES.items():
                res = run_dfcp(
                    args.seq_file,

                    val=0.2, mask=mask,

                    tree=args.tree, variant_pos_fname=args.variant_pos_fname,
                    variant_start_pos=args.variant_start_pos,
                    clade_beta=args.clade_beta,

                    block_init=int(init_method == "block"),
                    pbwt_init=int(init_method[0] == "pbwt"),
                    pbwt_match_len=init_method[1] if init_method[0] == "pbwt" else None,
                    init_only=init_only,
                )
                for metric in SHARED_METRICS:
                    data[phase][init_method][metric].append(res[metric])
                if phase == "post":
                    data[phase][init_method]["elbo"].append(res["train_log"][-1]["elbo"])
                    data[phase][init_method]["iterations"].append(len(res["train_log"]))

    ncols = 2
    nrows = (len(METRICS) + ncols - 1) // ncols
    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=tuple(PLOT_TITLES.get(metric, metric) for metric in METRICS),
        vertical_spacing=0.16, horizontal_spacing=0.08,
    )
    for plot_idx, metric in enumerate(METRICS):
        row = plot_idx // ncols + 1
        col = plot_idx % ncols + 1
        phases = ("post",) if metric in POST_METRICS else PHASES
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
                        y=data[phase][init_method][metric],
                        name=name,
                        legendgroup=name,
                        showlegend=plot_idx == 0,
                        visible=trace_phase in ("pre", "post_only"),
                        line=line,
                        meta={
                            "method": method,
                            "match_len": init_method[1] if is_pbwt else None,
                            "phase": trace_phase,
                        },
                    ),
                    row=row,
                    col=col,
                )

    sliders = []
    slider_steps = [
        {"label": str(match_len), "method": "skip"}
        for match_len in PBWT_MATCH_LENGTHS
    ]
    for plot_idx in range(len(METRICS)):
        axis_suffix = "" if plot_idx == 0 else str(plot_idx + 1)
        x_domain = fig.layout[f"xaxis{axis_suffix}"].domain
        y_domain = fig.layout[f"yaxis{axis_suffix}"].domain
        sliders.append({
            "active": 0,
            "currentvalue": {"prefix": "PBWT match length: "},
            "steps": slider_steps,
            "x": x_domain[0],
            "len": x_domain[1] - x_domain[0],
            "y": y_domain[0] - 0.04,
            "yanchor": "top",
        })

    fig.update_xaxes(title_text="mask frac")
    fig.update_layout(
        height=nrows * 500,
        title=f"DFCP init methods: pre/post-training metrics<br>seq_file: {args.seq_file}",
        margin={"b": 140, "t": 130},
        sliders=sliders,
        updatemenus=[
            {
                "active": 0,
                "buttons": [
                    {"label": "All PBWT", "method": "skip"},
                    {"label": "Selected PBWT", "method": "skip"},
                ],
                "direction": "left",
                "showactive": True,
                "type": "buttons",
                "x": 1,
                "xanchor": "right",
                "y": 1.06,
                "yanchor": "bottom",
            },
            {
                "active": 0,
                "buttons": [
                    {"label": "Pre-training", "method": "skip"},
                    {"label": "Post-training", "method": "skip"},
                ],
                "direction": "left",
                "showactive": True,
                "type": "buttons",
                "x": 0.65,
                "xanchor": "right",
                "y": 1.06,
                "yanchor": "bottom",
            },
        ],
    )

    post_script = r"""
const plot = document.getElementById('{plot_id}');
let selectedMatchLength = Number(plot.layout.sliders[0].steps[0].label);
let showAllPbwt = true;
let trainingPhase = 'pre';

function updateVisibleTraces() {
    const visible = plot.data.map(trace =>
        (trace.meta.phase === trainingPhase || trace.meta.phase === 'post_only') &&
        (trace.meta.method !== 'pbwt' || showAllPbwt ||
            trace.meta.match_len === selectedMatchLength)
    );
    Plotly.restyle(plot, {visible: visible});
}

plot.on('plotly_sliderchange', event => {
    selectedMatchLength = Number(event.step.label);
    showAllPbwt = false;
    const active = plot.layout.sliders[0].steps.findIndex(
        step => Number(step.label) === selectedMatchLength
    );
    const update = {'updatemenus[0].active': 1};
    plot.layout.sliders.forEach((slider, index) => {
        update[`sliders[${index}].active`] = active;
    });
    Plotly.relayout(plot, update);
    updateVisibleTraces();
});

plot.on('plotly_buttonclicked', event => {
    if (event.button.label === 'All PBWT' || event.button.label === 'Selected PBWT') {
        showAllPbwt = event.button.label === 'All PBWT';
    } else {
        trainingPhase = event.button.label === 'Pre-training' ? 'pre' : 'post';
    }
    updateVisibleTraces();
});
"""
    fig.write_html("init.html", post_script=post_script)
