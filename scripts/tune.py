import argparse
import itertools
import json
import math
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import plotly.graph_objects as go
import plotly.io as pio
import torch
from botorch.acquisition.logei import qLogNoisyExpectedImprovement
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms.outcome import Standardize
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from gpytorch.mlls import ExactMarginalLogLikelihood
from plotly.subplots import make_subplots


DEFAULT_SEQ_DIR = Path("data/examples/simulated/SIM1_LEN500_NHAPS100")
PBWT_MATCH_LEN = 20
CONFIRM_REPLICATES = 5
NEAR_BEST_TOLERANCE = 0.0025
RECOMMENDATION_POINTS = 8192
MC_SAMPLES = 128
ACQ_RESTARTS = 10
ACQ_RAW_SAMPLES = 256
FAILURE_RETRIES = 1
GRID_SIZE = 45
MODE_COLORS = {"hard": "#3366CC", "noisy": "#C9342C", "soft": "#2A9D8F"}


@dataclass(frozen=True)
class ParameterSpec:
    name: str
    label: str
    lower: float
    upper: float
    scale: str
    group: str

    def decode(self, unit_value: float) -> float:
        if self.scale == "log":
            return math.exp(
                math.log(self.lower)
                + unit_value * (math.log(self.upper) - math.log(self.lower))
            )
        return self.lower + unit_value * (self.upper - self.lower)

    def encode(self, value: float) -> float:
        if self.scale == "log":
            return (math.log(value) - math.log(self.lower)) / (
                math.log(self.upper) - math.log(self.lower)
            )
        return (value - self.lower) / (self.upper - self.lower)


COMMON_SPECS = (
    ParameterSpec("alpha_mean", "alpha prior mean", 0.05, 20.0, "log", "alpha"),
    ParameterSpec("alpha_strength", "alpha prior shape", 1.0, 8.0, "log", "alpha"),
    ParameterSpec("v_1", "discount prior v1", 1.0, 8.0, "log", "discount"),
    ParameterSpec("v_2", "discount prior v2", 1.0, 8.0, "log", "discount"),
    ParameterSpec("gamma_mean", "gamma prior mean", 0.01, 10.0, "log", "gamma"),
    ParameterSpec("gamma_strength", "gamma prior shape", 1.0, 8.0, "log", "gamma"),
)
NOISY_SPECS = (
    ParameterSpec("error_mean", "mismatch prior mean", 0.001, 0.3, "log", "error"),
    ParameterSpec(
        "error_strength", "mismatch prior concentration", 1.0, 10.0, "log", "error"
    ),
)


def get_specs(mode: str) -> tuple[ParameterSpec, ...]:
    return COMMON_SPECS + (NOISY_SPECS if mode == "noisy" else ())


def decode_point(x: list[float], specs: tuple[ParameterSpec, ...]) -> dict[str, float]:
    params = {spec.name: spec.decode(value) for spec, value in zip(specs, x)}
    params["discount_mean"] = params["v_1"] / (params["v_1"] + params["v_2"])
    params["discount_strength"] = params["v_1"] + params["v_2"]
    return params


def to_cpp_params(search_params: dict[str, float], mode: str) -> dict[str, float]:
    cpp_params = {
        "tau_1": search_params["alpha_strength"],
        "tau_2": search_params["alpha_strength"] / search_params["alpha_mean"],
        "v_1": search_params["v_1"],
        "v_2": search_params["v_2"],
        "phi_1": search_params["gamma_strength"],
        "phi_2": search_params["gamma_strength"] / search_params["gamma_mean"],
    }
    if mode == "noisy":
        cpp_params["lambda_1"] = (
            search_params["error_mean"] * search_params["error_strength"]
        )
        cpp_params["lambda_2"] = (
            (1.0 - search_params["error_mean"]) * search_params["error_strength"]
        )
    return cpp_params


def current_default_x(mode: str, specs: tuple[ParameterSpec, ...]) -> list[float]:
    values = {
        "alpha_mean": 1.0,
        "alpha_strength": 1.0,
        "v_1": 1.0,
        "v_2": 1.0,
        "gamma_mean": 1.0,
        "gamma_strength": 2.0,
        "error_mean": 0.1 / 1.1,
        "error_strength": 1.1,
    }
    return [spec.encode(values[spec.name]) for spec in specs]


def resolve_seq_files(args: argparse.Namespace) -> list[Path]:
    if args.seq_files:
        files = args.seq_files
    elif args.seq_dir:
        files = sorted(args.seq_dir.glob("haps*.txt*"))
    else:
        files = sorted(DEFAULT_SEQ_DIR.glob("haps*.txt*"))
    if not files:
        raise ValueError("no haplotype files were found")
    missing = [path for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing sequence file: {missing[0]}")
    return files


def build_dfcp() -> None:
    subprocess.run(["./build.sh"], check=True)


def get_dfcp_command(
    seq_file: Path, mode: str, cpp_params: dict[str, float], mask: float, val: float
) -> list[str]:
    command = [
        "./build/dfcp",
        str(seq_file),
        "--pbwt_init",
        "1",
        "--pbwt_match_len",
        str(PBWT_MATCH_LEN),
        "--mask",
        str(mask),
        "--val",
        str(val),
        "--soft",
        "1" if mode == "soft" else "0",
        "--noisy",
        "1" if mode == "noisy" else "0",
    ]
    for name, value in cpp_params.items():
        command.extend((f"--{name}", repr(value)))
    return command


def run_dfcp(command: list[str]) -> dict:
    last_error = ""
    for _ in range(FAILURE_RETRIES + 1):
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            try:
                return json.loads(result.stdout)
            except json.JSONDecodeError as error:
                last_error = f"invalid JSON: {error}"
        else:
            last_error = result.stderr.strip()[-1000:]
    raise RuntimeError(last_error or "DFCP failed without an error message")


def evaluate_candidate(
    x: list[float],
    specs: tuple[ParameterSpec, ...],
    mode: str,
    seq_files: list[Path],
    mask: float,
    val: float,
    replicates: int,
    source: str,
    iteration: int,
) -> dict:
    search_params = decode_point(x, specs)
    cpp_params = to_cpp_params(search_params, mode)
    replicate_scores = []
    within_replicate_variances = []
    runs = []
    for replicate in range(replicates):
        dataset_scores = []
        dataset_variances = []
        for seq_file in seq_files:
            output = run_dfcp(get_dfcp_command(seq_file, mode, cpp_params, mask, val))
            accuracy = float(output["dfcp_impute_acc"])
            n_masked = int(output["n_masked_alleles"])
            dataset_scores.append(accuracy)
            dataset_variances.append(accuracy * (1.0 - accuracy) / n_masked)
            runs.append({
                "replicate": replicate,
                "seq_file": str(seq_file),
                "accuracy": accuracy,
                "n_masked_alleles": n_masked,
            })
        replicate_scores.append(sum(dataset_scores) / len(dataset_scores))
        within_replicate_variances.append(
            sum(dataset_variances) / (len(dataset_variances) ** 2)
        )

    objective_mean = sum(replicate_scores) / replicates
    sample_variance = sum(
        (score - objective_mean) ** 2 for score in replicate_scores
    ) / (replicates - 1)
    variance_of_mean = sample_variance / replicates
    binomial_floor = sum(within_replicate_variances) / (replicates**2)
    objective_variance = max(variance_of_mean, binomial_floor, 1e-8)
    record = {
        "x": x,
        "search_params": search_params,
        "cpp_params": cpp_params,
        "objective_mean": objective_mean,
        "objective_variance": objective_variance,
        "replicate_scores": replicate_scores,
        "runs": runs,
        "source": source,
        "iteration": iteration,
    }
    print(
        f"{mode:5s} {source:28s} {iteration:3d} "
        f"accuracy={objective_mean:.6f} se={math.sqrt(objective_variance):.6f}",
        flush=True,
    )
    return record


def fit_surrogate(records: list[dict]) -> SingleTaskGP:
    train_x = torch.tensor([record["x"] for record in records], dtype=torch.double)
    train_y = torch.tensor(
        [[record["objective_mean"]] for record in records], dtype=torch.double
    )
    train_yvar = torch.tensor(
        [[record["objective_variance"]] for record in records], dtype=torch.double
    )
    model = SingleTaskGP(
        train_x,
        train_y,
        train_Yvar=train_yvar,
        outcome_transform=Standardize(m=1),
    )
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))
    return model


def next_bo_point(records: list[dict]) -> list[float]:
    model = fit_surrogate(records)
    train_x = torch.tensor([record["x"] for record in records], dtype=torch.double)
    acquisition = qLogNoisyExpectedImprovement(
        model=model,
        X_baseline=train_x,
        sampler=SobolQMCNormalSampler(sample_shape=torch.Size([MC_SAMPLES])),
    )
    bounds = torch.stack(
        (
            torch.zeros(train_x.shape[1], dtype=torch.double),
            torch.ones(train_x.shape[1], dtype=torch.double),
        )
    )
    candidate, _ = optimize_acqf(
        acquisition,
        bounds=bounds,
        q=1,
        num_restarts=ACQ_RESTARTS,
        raw_samples=ACQ_RAW_SAMPLES,
    )
    return candidate.squeeze(0).tolist()


def less_informative_penalty(x: torch.Tensor, specs: tuple[ParameterSpec, ...]) -> float:
    params = decode_point(x.tolist(), specs)
    penalty = (
        math.log(params["alpha_strength"])
        + math.log(params["gamma_strength"])
        + abs(math.log(params["v_1"]))
        + abs(math.log(params["v_2"]))
    )
    if "error_strength" in params:
        penalty += math.log(params["error_strength"])
    return penalty


def select_candidates(
    model: SingleTaskGP, specs: tuple[ParameterSpec, ...], mode: str, seed: int
) -> tuple[dict, dict]:
    candidates = torch.quasirandom.SobolEngine(
        len(specs), scramble=True, seed=seed
    ).draw(RECOMMENDATION_POINTS).double()
    with torch.no_grad():
        posterior_means = model.posterior(candidates).mean.squeeze(-1)
    best_index = int(torch.argmax(posterior_means))
    best_mean = float(posterior_means[best_index])
    eligible = torch.nonzero(
        posterior_means >= best_mean - NEAR_BEST_TOLERANCE
    ).squeeze(-1)
    less_informative_index = min(
        eligible.tolist(), key=lambda index: less_informative_penalty(candidates[index], specs)
    )

    def describe(index: int) -> dict:
        x = candidates[index].tolist()
        search_params = decode_point(x, specs)
        return {
            "x": x,
            "posterior_mean": float(posterior_means[index]),
            "search_params": search_params,
            "cpp_params": to_cpp_params(search_params, mode),
        }

    return describe(best_index), describe(less_informative_index)


def get_ard_summary(model: SingleTaskGP, specs: tuple[ParameterSpec, ...]) -> tuple[dict, dict]:
    lengthscales = model.covar_module.lengthscale.detach().reshape(-1).tolist()
    raw = {spec.name: value for spec, value in zip(specs, lengthscales)}
    inverse = {name: 1.0 / value for name, value in raw.items()}
    total = sum(inverse.values())
    return raw, {name: value / total for name, value in inverse.items()}


def tune_mode(
    mode: str, seq_files: list[Path], args: argparse.Namespace
) -> dict:
    specs = get_specs(mode)
    default_record = evaluate_candidate(
        current_default_x(mode, specs), specs, mode, seq_files, args.mask, args.val,
        CONFIRM_REPLICATES, "current_default", -1,
    )
    records = []
    sobol = torch.quasirandom.SobolEngine(
        len(specs), scramble=True, seed=args.seed + 1000 * args.modes.index(mode)
    )
    for iteration, point in enumerate(sobol.draw(args.initial_points).double()):
        records.append(evaluate_candidate(
            point.tolist(), specs, mode, seq_files, args.mask, args.val,
            args.replicates, "sobol", iteration,
        ))
    for iteration in range(args.iterations):
        point = next_bo_point(records)
        records.append(evaluate_candidate(
            point, specs, mode, seq_files, args.mask, args.val,
            args.replicates, "bo", iteration,
        ))

    selection_model = fit_surrogate(records)
    predicted_optimum, less_informative = select_candidates(
        selection_model, specs, mode, args.seed + 10000 + args.modes.index(mode)
    )

    predicted_confirmation = evaluate_candidate(
        predicted_optimum["x"], specs, mode, seq_files, args.mask, args.val,
        CONFIRM_REPLICATES, "confirm_predicted_optimum", len(records),
    )
    less_informative_confirmation = evaluate_candidate(
        less_informative["x"], specs, mode, seq_files, args.mask, args.val,
        CONFIRM_REPLICATES, "confirm_less_informative", len(records) + 1,
    )
    predicted_optimum["confirmation"] = predicted_confirmation
    less_informative["confirmation"] = less_informative_confirmation
    records.extend((predicted_confirmation, less_informative_confirmation))

    final_model = fit_surrogate(records)
    lengthscales, normalized_inverse = get_ard_summary(final_model, specs)
    core_records = [record for record in records if record["source"] in ("sobol", "bo")]
    return {
        "parameter_specs": [asdict(spec) for spec in specs],
        "records": records,
        "summary": {
            "current_default": default_record,
            "best_observed": max(core_records, key=lambda record: record["objective_mean"]),
            "predicted_optimum": predicted_optimum,
            "less_informative": less_informative,
            "ard_lengthscales": lengthscales,
            "normalized_inverse_lengthscales": normalized_inverse,
            "lengthscale_caution": (
                "ARD lengthscales describe this fitted surrogate and sampled region; "
                "they are not causal importance measures."
            ),
        },
    }


def get_pairs(specs: tuple[ParameterSpec, ...]) -> list[tuple[int, int]]:
    pairs = list(itertools.combinations(range(len(specs)), 2))
    return sorted(pairs, key=lambda pair: (specs[pair[0]].group != specs[pair[1]].group, pair))


def predict_in_chunks(
    model: SingleTaskGP, points: torch.Tensor, chunk_size: int = 4096
) -> tuple[torch.Tensor, torch.Tensor]:
    means = []
    variances = []
    with torch.no_grad():
        for chunk in points.split(chunk_size):
            posterior = model.posterior(chunk)
            means.append(posterior.mean.squeeze(-1))
            variances.append(posterior.variance.squeeze(-1))
    return torch.cat(means), torch.cat(variances)


def get_surface(
    model: SingleTaskGP,
    anchor: torch.Tensor,
    specs: tuple[ParameterSpec, ...],
    pair: tuple[int, int],
) -> dict:
    unit_axis = torch.linspace(0.0, 1.0, GRID_SIZE, dtype=torch.double)
    x_unit, y_unit = torch.meshgrid(unit_axis, unit_axis, indexing="xy")
    points = anchor.repeat(GRID_SIZE * GRID_SIZE, 1)
    points[:, pair[0]] = x_unit.reshape(-1)
    points[:, pair[1]] = y_unit.reshape(-1)
    mean, variance = predict_in_chunks(model, points)
    return {
        "x": [specs[pair[0]].decode(float(value)) for value in unit_axis],
        "y": [specs[pair[1]].decode(float(value)) for value in unit_axis],
        "z": mean.reshape(GRID_SIZE, GRID_SIZE).tolist(),
        "std": variance.clamp_min(0.0).sqrt().reshape(GRID_SIZE, GRID_SIZE).tolist(),
    }


def get_mode_report(mode: str, mode_data: dict) -> dict:
    specs = tuple(ParameterSpec(**spec) for spec in mode_data["parameter_specs"])
    records = mode_data["records"]
    model = fit_surrogate(records)
    anchor = torch.tensor(mode_data["summary"]["less_informative"]["x"], dtype=torch.double)
    pair_reports = []
    for pair in get_pairs(specs):
        pair_reports.append({
            "id": f"{specs[pair[0]].name}|{specs[pair[1]].name}",
            "x_name": specs[pair[0]].name,
            "x_label": specs[pair[0]].label,
            "y_name": specs[pair[1]].name,
            "y_label": specs[pair[1]].label,
            "surface": get_surface(model, anchor, specs, pair),
            "observed_x": [record["search_params"][specs[pair[0]].name] for record in records],
            "observed_y": [record["search_params"][specs[pair[1]].name] for record in records],
            "observed_accuracy": [record["objective_mean"] for record in records],
        })
    surface_values = [
        value for pair in pair_reports for row in pair["surface"]["z"] for value in row
    ]
    accuracies = [record["objective_mean"] for record in records]
    return {
        "mode": mode,
        "specs": specs,
        "records": records,
        "pairs": pair_reports,
        "zmin": min(surface_values + accuracies),
        "zmax": max(surface_values + accuracies),
        "summary": mode_data["summary"],
    }


def add_mode_traces(fig: go.Figure, report: dict) -> None:
    mode = report["mode"]
    sensitivity = report["summary"]["normalized_inverse_lengthscales"]
    ordered = sorted(sensitivity, key=sensitivity.get)
    labels = {spec.name: spec.label for spec in report["specs"]}
    fig.add_trace(go.Bar(
        x=[sensitivity[name] for name in ordered], y=[labels[name] for name in ordered],
        orientation="h", marker_color=MODE_COLORS[mode], showlegend=False,
        hovertemplate="%{y}<br>normalized inverse lengthscale=%{x:.3f}<extra></extra>",
        meta={"kind": "sensitivity"},
    ), row=1, col=1)

    first_pair = report["pairs"][0]["id"]
    for pair in report["pairs"]:
        visible = pair["id"] == first_pair
        meta = {"kind": "surface", "pair": pair["id"]}
        fig.add_trace(go.Heatmap(
            x=pair["surface"]["x"], y=pair["surface"]["y"], z=pair["surface"]["z"],
            customdata=pair["surface"]["std"], coloraxis="coloraxis",
            hovertemplate=(
                f"{pair['x_label']}: %{{x:.4g}}<br>{pair['y_label']}: %{{y:.4g}}<br>"
                "predicted accuracy: %{z:.5f}<br>posterior SD: %{customdata:.5f}<extra></extra>"
            ), showlegend=False, meta=meta, visible=visible,
        ), row=1, col=2)
        fig.add_trace(go.Contour(
            x=pair["surface"]["x"], y=pair["surface"]["y"], z=pair["surface"]["z"],
            contours={"coloring": "none", "showlabels": True},
            line={"color": "rgba(255,255,255,0.8)", "width": 1},
            showlegend=False, showscale=False, hoverinfo="skip", meta=meta,
            visible=visible,
        ), row=1, col=2)
        fig.add_trace(go.Scattergl(
            x=pair["observed_x"], y=pair["observed_y"], mode="markers",
            marker={"color": pair["observed_accuracy"], "coloraxis": "coloraxis",
                    "line": {"color": "white", "width": 1}, "size": 7},
            text=[f"accuracy={value:.5f}" for value in pair["observed_accuracy"]],
            hovertemplate="%{text}<extra>evaluated point</extra>", showlegend=False,
            meta=meta, visible=visible,
        ), row=1, col=2)
        optimum = report["summary"]["predicted_optimum"]
        fig.add_trace(go.Scatter(
            x=[optimum["search_params"][pair["x_name"]]],
            y=[optimum["search_params"][pair["y_name"]]],
            mode="markers",
            marker={
                "symbol": "star",
                "size": 16,
                "color": "black",
                "line": {"color": "white", "width": 1.5},
            },
            hovertemplate=(
                "GP optimum candidate<br>"
                f"{pair['x_label']}: %{{x:.4g}}<br>"
                f"{pair['y_label']}: %{{y:.4g}}<extra></extra>"
            ),
            showlegend=False,
            meta=meta,
            visible=visible,
        ), row=1, col=2)

    core = [record for record in report["records"] if record["source"] in ("sobol", "bo")]
    running_best = -math.inf
    best = []
    for record in core:
        running_best = max(running_best, record["objective_mean"])
        best.append(running_best)
    fig.add_trace(go.Scatter(
        x=list(range(1, len(core) + 1)), y=best, mode="lines+markers",
        line={"color": MODE_COLORS[mode]}, showlegend=False,
        meta={"kind": "convergence"},
        hovertemplate="evaluation %{x}<br>best accuracy=%{y:.5f}<extra></extra>",
    ), row=2, col=1)
    fig.add_trace(go.Scatter(
        x=list(range(1, len(report["records"]) + 1)),
        y=[record["objective_mean"] for record in report["records"]],
        error_y={"type": "data", "array": [
            math.sqrt(record["objective_variance"]) for record in report["records"]
        ], "visible": True},
        text=[record["source"] for record in report["records"]], mode="markers",
        marker={"color": MODE_COLORS[mode], "size": 7}, showlegend=False,
        meta={"kind": "observations"},
        hovertemplate="evaluation %{x}<br>accuracy=%{y:.5f}<br>%{text}<extra></extra>",
    ), row=2, col=2)


def get_sidebar_context(reports: list[dict]) -> dict:
    context = {}
    for report in reports:
        summary = report["summary"]
        context[report["mode"]] = {
            "pairs": [{
                "id": pair["id"],
                "label": f"{pair['x_label']} × {pair['y_label']}",
                "x_label": pair["x_label"],
                "y_label": pair["y_label"],
                "x_scale": next(
                    spec.scale for spec in report["specs"] if spec.name == pair["x_name"]
                ),
                "y_scale": next(
                    spec.scale for spec in report["specs"] if spec.name == pair["y_name"]
                ),
            } for pair in report["pairs"]],
            "best_accuracy": summary["best_observed"]["objective_mean"],
            "best_runs": len(summary["best_observed"]["runs"]),
            "current_default_accuracy": summary["current_default"]["objective_mean"],
            "current_default_runs": len(summary["current_default"]["runs"]),
            "predicted_optimum": summary["predicted_optimum"],
            "less_informative": summary["less_informative"],
        }
    return context


def get_page_script(context: dict) -> str:
    template = r"""
const report = __REPORT_CONTEXT__;
let selectedMode = Object.keys(report)[0];
let selectedPair = report[selectedMode].pairs[0].id;

function updateSummary() {
    const mode = report[selectedMode];
    document.getElementById('current-default-label').textContent =
        `Current defaults: actual mean (${mode.current_default_runs} runs)`;
    document.getElementById('current-default-accuracy').textContent =
        mode.current_default_accuracy.toFixed(5);
    document.getElementById('best-label').textContent =
        `Best observed: actual mean (${mode.best_runs} runs)`;
    document.getElementById('best-accuracy').textContent = mode.best_accuracy.toFixed(5);
    document.getElementById('predicted-accuracy').textContent =
        mode.predicted_optimum.posterior_mean.toFixed(5);
    document.getElementById('predicted-actual-label').textContent =
        `Optimum: actual mean (${mode.predicted_optimum.confirmation.runs.length} runs)`;
    document.getElementById('predicted-confirmed-accuracy').textContent =
        mode.predicted_optimum.confirmation.objective_mean.toFixed(5);
    document.getElementById('less-informative-accuracy').textContent =
        mode.less_informative.posterior_mean.toFixed(5);
    document.getElementById('less-informative-actual-label').textContent =
        `Less-informative: actual mean (${mode.less_informative.confirmation.runs.length} runs)`;
    document.getElementById('less-informative-confirmed-accuracy').textContent =
        mode.less_informative.confirmation.objective_mean.toFixed(5);
    document.getElementById('candidate-params').innerHTML = Object.entries(
        mode.less_informative.cpp_params
    ).map(([name, value]) =>
        `<tr><td>${name}</td><td>${Number(value).toPrecision(4)}</td></tr>`
    ).join('');
    document.getElementById('optimum-params').innerHTML = Object.entries(
        mode.predicted_optimum.cpp_params
    ).map(([name, value]) =>
        `<tr><td>${name}</td><td>${Number(value).toPrecision(4)}</td></tr>`
    ).join('');
}

function updatePairOptions() {
    const select = document.getElementById('pair-select');
    select.innerHTML = report[selectedMode].pairs.map(pair =>
        `<option value="${pair.id}">${pair.label}</option>`
    ).join('');
    selectedPair = report[selectedMode].pairs[0].id;
    select.value = selectedPair;
}

function updatePlot() {
    document.querySelectorAll('.mode-report').forEach(section => {
        section.hidden = section.dataset.mode !== selectedMode;
    });
    const plot = document.getElementById(`tune-plot-${selectedMode}`);
    const visible = plot.data.map(trace =>
        trace.meta.kind !== 'surface' || trace.meta.pair === selectedPair
    );
    const pair = report[selectedMode].pairs.find(candidate => candidate.id === selectedPair);
    Plotly.restyle(plot, {visible: visible});
    Plotly.relayout(plot, {
        'xaxis2.title.text': pair.x_label,
        'yaxis2.title.text': pair.y_label,
        'xaxis2.type': pair.x_scale,
        'yaxis2.type': pair.y_scale,
        'xaxis2.autorange': true,
        'yaxis2.autorange': true,
    });
    updateSummary();
    requestAnimationFrame(() => Plotly.Plots.resize(plot));
}

const sidebar = document.createElement('aside');
sidebar.id = 'tune-sidebar';
sidebar.innerHTML = `
    <h2>Hyperparameter report</h2>
    <section><h3>DFCP mode</h3><div class="tune-buttons" id="mode-controls">
        ${Object.keys(report).map((mode, index) =>
            `<button class="${index === 0 ? 'active' : ''}" data-mode="${mode}">${mode}</button>`
        ).join('')}
    </div></section>
    <section><h3>Parameter pair</h3><select id="pair-select"></select></section>
    <section class="tune-summary"><h3>Accuracy summary</h3>
        <div><span id="current-default-label"></span><strong id="current-default-accuracy"></strong></div>
        <div><span id="best-label"></span><strong id="best-accuracy"></strong></div>
        <div><span>Optimum: GP prediction</span><strong id="predicted-accuracy"></strong></div>
        <div><span id="predicted-actual-label"></span><strong id="predicted-confirmed-accuracy"></strong></div>
        <div><span>Less-informative: GP prediction</span><strong id="less-informative-accuracy"></strong></div>
        <div><span id="less-informative-actual-label"></span><strong id="less-informative-confirmed-accuracy"></strong></div>
    </section>
    <section><h3>GP optimum C++ candidate</h3>
        <table><tbody id="optimum-params"></tbody></table></section>
    <section><h3>Less-informative C++ candidate</h3>
        <table><tbody id="candidate-params"></tbody></table></section>
    <section class="tune-note"><h3>Interpretation</h3>
        <p>The heatmap is the GP posterior mean with other parameters fixed at the
        less-informative candidate. White contours show equal predicted accuracy.</p>
        <p>The black star is the GP optimum. Dots are evaluated configurations
        projected onto the selected pair; their other parameters varied.</p>
        <p>Shorter ARD lengthscales indicate faster surrogate variation, but are not
        causal importance scores.</p>
    </section>`;
document.body.prepend(sidebar);

const style = document.createElement('style');
style.textContent = `
    html, body { margin: 0; }
    #tune-sidebar { position: fixed; inset: 0 auto 0 0; z-index: 1000;
        box-sizing: border-box; width: 350px; padding: 18px; overflow-y: auto;
        border-right: 1px solid #d9dce1; background: rgba(255,255,255,.98);
        color: #2a3f5f; font-family: Arial,sans-serif;
        box-shadow: 2px 0 8px rgba(0,0,0,.08); }
    #tune-sidebar h2 { margin: 0 0 18px; font-size: 20px; }
    #tune-sidebar h3 { margin: 0 0 8px; font-size: 14px; }
    #tune-sidebar section { margin-bottom: 20px; }
    #pair-select { box-sizing: border-box; width: 100%; padding: 7px 5px;
        border: 1px solid #aeb5c0; border-radius: 4px; background: white; color: #2a3f5f; }
    .tune-buttons { display: flex; gap: 5px; }
    .tune-buttons button { flex: 1; padding: 7px 5px; border: 1px solid #aeb5c0;
        border-radius: 4px; background: #f4f5f7; color: #2a3f5f; cursor: pointer;
        text-transform: capitalize; }
    .tune-buttons button.active { border-color: #3366cc; background: #3366cc; color: white; }
    .tune-summary > div { display: grid; grid-template-columns: minmax(0,1fr) auto;
        gap: 10px; margin: 8px 0; font-size: 12px; line-height: 1.25; }
    #tune-sidebar table { width: 100%; border-collapse: collapse; font-size: 12px; }
    #tune-sidebar td { padding: 4px 2px; border-bottom: 1px solid #e5e7eb; }
    #tune-sidebar td:last-child { text-align: right; font-family: monospace; }
    .tune-note { color: #5d6778; font-size: 12px; line-height: 1.4; }
    .mode-report { width: calc(100% - 350px); margin-left: 350px; }
    .mode-report[hidden] { display: none; }
    @media (max-width: 850px) { #tune-sidebar { width: 280px; padding: 12px; }
        .mode-report { width: calc(100% - 280px); margin-left: 280px; } }`;
document.head.appendChild(style);

document.querySelectorAll('#mode-controls button').forEach(button => {
    button.addEventListener('click', () => {
        document.querySelectorAll('#mode-controls button').forEach(
            other => other.classList.remove('active'));
        button.classList.add('active');
        selectedMode = button.dataset.mode;
        updatePairOptions();
        updatePlot();
    });
});
document.getElementById('pair-select').addEventListener('change', event => {
    selectedPair = event.target.value;
    updatePlot();
});
updatePairOptions();
updatePlot();
window.addEventListener('resize', () => {
    Plotly.Plots.resize(document.getElementById(`tune-plot-${selectedMode}`));
});
"""
    return template.replace("__REPORT_CONTEXT__", json.dumps(context))


def make_mode_figure(report: dict) -> go.Figure:
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "ARD sensitivity (normalized inverse lengthscale)",
            "Pairwise GP posterior accuracy",
            "Best observed accuracy",
            "Observed accuracy and estimated noise",
        ),
        vertical_spacing=0.16,
        horizontal_spacing=0.15,
    )
    add_mode_traces(fig, report)
    first_pair = report["pairs"][0]
    fig.update_xaxes(title_text=first_pair["x_label"], type="log", row=1, col=2)
    fig.update_yaxes(title_text=first_pair["y_label"], type="log", row=1, col=2)
    fig.update_xaxes(title_text="BO evaluation", row=2, col=1)
    fig.update_xaxes(title_text="evaluation", row=2, col=2)
    fig.update_yaxes(title_text="accuracy", row=2, col=1)
    fig.update_yaxes(title_text="accuracy", row=2, col=2)
    fig.update_layout(
        height=980,
        margin={"t": 115, "b": 70, "l": 95, "r": 90},
        title=(
            "DFCP hyperparameter optimization and sensitivity"
            f" — {report['mode']} mode"
            f"<br>PBWT initialization, match length {PBWT_MATCH_LEN}"
        ),
        showlegend=False,
        coloraxis={
            "colorscale": "Viridis",
            "cmin": report["zmin"],
            "cmax": report["zmax"],
            "colorbar": {
                "title": "predicted accuracy",
                "x": 1.015,
                "xanchor": "left",
                "y": 0.79,
                "yanchor": "middle",
                "len": 0.42,
            },
        },
    )
    return fig


def write_report(experiment: dict, output: Path) -> None:
    reports = [
        get_mode_report(mode, mode_data) for mode, mode_data in experiment["modes"].items()
    ]
    fragments = []
    for index, report in enumerate(reports):
        figure = make_mode_figure(report)
        plot = pio.to_html(
            figure,
            full_html=False,
            include_plotlyjs=True if index == 0 else False,
            config={"responsive": True},
            div_id=f"tune-plot-{report['mode']}",
        )
        hidden = "" if index == 0 else " hidden"
        fragments.append(
            f'<section class="mode-report" data-mode="{report["mode"]}"{hidden}>'
            f"{plot}</section>"
        )

    page = (
        '<!doctype html>\n<html><head><meta charset="utf-8">'
        "<title>DFCP hyperparameter report</title></head><body>\n"
        + "\n".join(fragments)
        + "\n<script>\n"
        + get_page_script(get_sidebar_context(reports))
        + "\n</script>\n</body></html>\n"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Tune hard, noisy, and soft DFCP priors with BoTorch and write a Plotly report."
        )
    )
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--seq_files", type=Path, nargs="+")
    inputs.add_argument(
        "--seq_dir", type=Path,
        help=(
            "use every haps*.txt* file in this directory; by default, use all such "
            f"files in {DEFAULT_SEQ_DIR}"
        ),
    )
    parser.add_argument(
        "--modes", nargs="+", choices=("hard", "noisy", "soft"),
        default=["hard", "noisy", "soft"],
    )
    parser.add_argument("--mask", type=float, default=0.2)
    parser.add_argument("--val", type=float, default=0.2)
    parser.add_argument("--initial_points", type=int, default=16)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=Path("output/tune2.html"))
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if len(set(args.modes)) != len(args.modes):
        raise ValueError("modes must not contain duplicates")
    if not 0.0 < args.mask <= 1.0 or not 0.0 < args.val < 1.0:
        raise ValueError("mask must be in (0,1] and val must be in (0,1)")
    if args.initial_points < 2:
        raise ValueError("initial_points must be at least 2")
    if args.iterations < 0:
        raise ValueError("iterations must be nonnegative")
    if args.replicates < 2:
        raise ValueError("replicates must be at least 2 to estimate observation noise")


def main() -> None:
    args = parse_args()
    validate_args(args)
    torch.manual_seed(args.seed)
    seq_files = resolve_seq_files(args)
    build_dfcp()
    experiment = {
        "version": 2,
        "config": {
            "seq_files": [str(path) for path in seq_files],
            "modes": args.modes,
            "mask": args.mask,
            "val": args.val,
            "pbwt_match_len": PBWT_MATCH_LEN,
            "initial_points": args.initial_points,
            "iterations": args.iterations,
            "replicates": args.replicates,
            "confirmation_replicates": CONFIRM_REPLICATES,
            "seed": args.seed,
            "objective": "dfcp_impute_acc",
            "parameterization": (
                "Gamma mean/shape, discount Beta v1/v2, mismatch Beta mean/concentration"
            ),
        },
        "modes": {},
    }
    for mode in args.modes:
        experiment["modes"][mode] = tune_mode(mode, seq_files, args)

    write_report(experiment, args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
