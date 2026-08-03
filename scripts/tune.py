import argparse
import itertools
import json
import math
import re
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

from plotly_html import ensure_plotly_asset


DEFAULT_SEQ_DIR = Path("data/examples/simulated/SIM1_LEN500_NHAPS100")
DEFAULT_SEQ_FILE = DEFAULT_SEQ_DIR / "haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt"
PBWT_MATCH_LEN = 20
CONFIRM_REPLICATES = 5
NEAR_BEST_TOLERANCE = 0.0025
RECOMMENDATION_POINTS = 8192
MC_SAMPLES = 128
ACQ_RESTARTS = 10
ACQ_RAW_SAMPLES = 256
FAILURE_RETRIES = 1
GRID_SIZE = 30
VISUAL_SIGNIFICANT_DIGITS = 7
PAIR_CELL_SIZE = 135
MODE_COLORS = {"hard": "#3366CC", "noisy": "#C9342C", "soft": "#2A9D8F"}
ERROR_RE = re.compile(r"\.txt\.gz_([0-9.]+)_([0-9.]+)\.txt\.gz_")


def get_seq_label(seq_file: Path) -> str:
    match = ERROR_RE.search(seq_file.name)
    if match is None:
        return "Baseline (no injected errors)"
    return f"Bit flip {match.group(1)}, switch {match.group(2)}"


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
        files = [DEFAULT_SEQ_FILE]
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
    seq_file: Path,
    mask: float,
    val: float,
    replicates: int,
    source: str,
    iteration: int,
) -> dict:
    search_params = decode_point(x, specs)
    cpp_params = to_cpp_params(search_params, mode)
    replicate_scores = []
    binomial_variances = []
    runs = []
    for replicate in range(replicates):
        output = run_dfcp(get_dfcp_command(seq_file, mode, cpp_params, mask, val))
        accuracy = float(output["dfcp_impute_acc"])
        n_masked = int(output["n_masked_alleles"])
        replicate_scores.append(accuracy)
        binomial_variances.append(accuracy * (1.0 - accuracy) / n_masked)
        runs.append({
            "replicate": replicate,
            "seq_file": str(seq_file),
            "accuracy": accuracy,
            "n_masked_alleles": n_masked,
        })

    objective_mean = sum(replicate_scores) / replicates
    sample_variance = sum(
        (score - objective_mean) ** 2 for score in replicate_scores
    ) / (replicates - 1)
    variance_of_mean = sample_variance / replicates
    binomial_floor = sum(binomial_variances) / (replicates**2)
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
        f"{seq_file.name} {mode:5s} {source:28s} {iteration:3d} "
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
    mode: str, seq_file: Path, dataset_index: int, args: argparse.Namespace
) -> dict:
    specs = get_specs(mode)
    mode_index = args.modes.index(mode)
    mode_seed = args.seed + 100_000 * dataset_index + 1_000 * mode_index
    default_record = evaluate_candidate(
        current_default_x(mode, specs), specs, mode, seq_file, args.mask, args.val,
        CONFIRM_REPLICATES, "current_default", -1,
    )
    records = []
    sobol = torch.quasirandom.SobolEngine(
        len(specs), scramble=True, seed=mode_seed
    )
    for iteration, point in enumerate(sobol.draw(args.initial_points).double()):
        records.append(evaluate_candidate(
            point.tolist(), specs, mode, seq_file, args.mask, args.val,
            args.replicates, "sobol", iteration,
        ))
    for iteration in range(args.iterations):
        point = next_bo_point(records)
        records.append(evaluate_candidate(
            point, specs, mode, seq_file, args.mask, args.val,
            args.replicates, "bo", iteration,
        ))

    selection_model = fit_surrogate(records)
    predicted_optimum, less_informative = select_candidates(
        selection_model, specs, mode, mode_seed + 10_000
    )

    predicted_confirmation = evaluate_candidate(
        predicted_optimum["x"], specs, mode, seq_file, args.mask, args.val,
        CONFIRM_REPLICATES, "confirm_predicted_optimum", len(records),
    )
    less_informative_confirmation = evaluate_candidate(
        less_informative["x"], specs, mode, seq_file, args.mask, args.val,
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


def round_visual(value: float) -> float:
    if value == 0.0:
        return 0.0
    places = VISUAL_SIGNIFICANT_DIGITS - 1 - math.floor(math.log10(abs(value)))
    return round(value, places)


def round_visual_grid(values: list[list[float]]) -> list[list[float]]:
    return [[round_visual(value) for value in row] for row in values]


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
        "x": [
            round_visual(specs[pair[0]].decode(float(value))) for value in unit_axis
        ],
        "y": [
            round_visual(specs[pair[1]].decode(float(value))) for value in unit_axis
        ],
        "z": round_visual_grid(mean.reshape(GRID_SIZE, GRID_SIZE).tolist()),
        "std": round_visual_grid(
            variance.clamp_min(0.0).sqrt().reshape(GRID_SIZE, GRID_SIZE).tolist()
        ),
    }


def get_mode_report(
    dataset_index: int, dataset: dict, mode: str, mode_data: dict
) -> dict:
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
            "observed_x": [
                round_visual(record["search_params"][specs[pair[0]].name])
                for record in records
            ],
            "observed_y": [
                round_visual(record["search_params"][specs[pair[1]].name])
                for record in records
            ],
            "observed_accuracy": [
                round_visual(record["objective_mean"]) for record in records
            ],
        })
    surface_values = [
        value for pair in pair_reports for row in pair["surface"]["z"] for value in row
    ]
    accuracies = [record["objective_mean"] for record in records]
    return {
        "dataset_id": f"d{dataset_index}",
        "dataset_label": dataset["label"],
        "seq_file": dataset["seq_file"],
        "mode": mode,
        "specs": specs,
        "records": records,
        "pairs": pair_reports,
        "zmin": min(surface_values + accuracies),
        "zmax": max(surface_values + accuracies),
        "summary": mode_data["summary"],
    }


def make_overview_figure(report: dict, width: int) -> go.Figure:
    mode = report["mode"]
    sensitivity = report["summary"]["normalized_inverse_lengthscales"]
    ordered = sorted(sensitivity, key=sensitivity.get)
    labels = {spec.name: spec.label for spec in report["specs"]}
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=(
            "ARD sensitivity (normalized inverse lengthscale)",
            "Selected pairwise GP posterior trend",
        ),
        horizontal_spacing=0.16,
    )
    fig.add_trace(go.Bar(
        x=[sensitivity[name] for name in ordered], y=[labels[name] for name in ordered],
        orientation="h", marker_color=MODE_COLORS[mode], showlegend=False,
        hovertemplate="%{y}<br>normalized inverse lengthscale=%{x:.3f}<extra></extra>",
        meta={"kind": "sensitivity"},
    ), row=1, col=1)

    fig.add_trace(go.Heatmap(
        x=[], y=[], z=[], coloraxis="coloraxis", showlegend=False,
        meta={"kind": "overview_surface", "role": "heatmap"},
    ), row=1, col=2)
    fig.add_trace(go.Contour(
        x=[], y=[], z=[], contours={"coloring": "none", "showlabels": True},
        showlegend=False, showscale=False,
        meta={"kind": "overview_surface", "role": "contour"},
    ), row=1, col=2)
    fig.add_trace(go.Scattergl(
        x=[], y=[], mode="markers", showlegend=False,
        meta={"kind": "overview_surface", "role": "observed"},
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=[], y=[], mode="markers", showlegend=False,
        meta={"kind": "overview_surface", "role": "optimum"},
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=[], y=[], mode="markers", showlegend=False,
        meta={"kind": "overview_surface", "role": "less_informative"},
    ), row=1, col=2)

    first_pair = report["pairs"][0]
    fig.update_xaxes(title_text="normalized inverse lengthscale", row=1, col=1)
    fig.update_xaxes(title_text=first_pair["x_label"], type="log", row=1, col=2)
    fig.update_yaxes(title_text=first_pair["y_label"], type="log", row=1, col=2)
    fig.update_layout(
        width=width,
        height=520,
        margin={"t": 115, "b": 70, "l": 180, "r": 105},
        showlegend=False,
        coloraxis={
            "colorscale": "Viridis",
            "cmin": report["zmin"],
            "cmax": report["zmax"],
            "colorbar": {"title": "predicted accuracy", "x": 1.02, "len": 0.7},
        },
    )
    return fig


def make_pair_grid_figure(report: dict, cell_size: int = PAIR_CELL_SIZE) -> go.Figure:
    specs = report["specs"]
    n_params = len(specs)
    grid_size = n_params - 1
    plot_size = cell_size * grid_size
    figure_size = plot_size + 250
    fig = make_subplots(
        rows=grid_size,
        cols=grid_size,
        horizontal_spacing=0.012,
        vertical_spacing=0.012,
    )
    for row in range(1, grid_size + 1):
        for col in range(1, grid_size + 1):
            fig.update_xaxes(visible=False, row=row, col=col)
            fig.update_yaxes(visible=False, row=row, col=col)

    param_index = {spec.name: index for index, spec in enumerate(specs)}
    optimum = report["summary"]["predicted_optimum"]
    less_informative = report["summary"]["less_informative"]
    def trace_meta(pair_id: str, role: str) -> dict[str, str]:
        return {"kind": "grid_surface", "pair": pair_id, "role": role}

    for pair in report["pairs"]:
        col = n_params - param_index[pair["y_name"]]
        row = param_index[pair["x_name"]] + 1
        z = [list(values) for values in zip(*pair["surface"]["z"])]
        std = [list(values) for values in zip(*pair["surface"]["std"])]
        fig.add_trace(go.Heatmap(
            x=pair["surface"]["y"], y=pair["surface"]["x"], z=z,
            customdata=std, coloraxis="coloraxis",
            hovertemplate=(
                f"{pair['y_label']}: %{{x:.4g}}<br>{pair['x_label']}: %{{y:.4g}}<br>"
                "predicted accuracy: %{z:.5f}<br>posterior SD: %{customdata:.5f}<extra></extra>"
            ), showlegend=False, meta=trace_meta(pair["id"], "heatmap"),
        ), row=row, col=col)
        fig.add_trace(go.Contour(
            x=pair["surface"]["y"], y=pair["surface"]["x"], z=z,
            contours={"coloring": "none", "showlabels": True},
            line={"color": "rgba(255,255,255,0.8)", "width": 1},
            textfont={"color": "white", "size": 9},
            showlegend=False, showscale=False, hoverinfo="skip",
            meta=trace_meta(pair["id"], "contour"),
        ), row=row, col=col)
        fig.add_trace(go.Scattergl(
            x=pair["observed_y"], y=pair["observed_x"], mode="markers",
            marker={"color": pair["observed_accuracy"], "coloraxis": "coloraxis",
                    "line": {"color": "black", "width": 0.5}, "size": 4,
                    "opacity": 1.0},
            text=[f"accuracy={value:.5f}" for value in pair["observed_accuracy"]],
            hovertemplate="%{text}<extra>evaluated point</extra>", showlegend=False,
            meta=trace_meta(pair["id"], "observed"),
        ), row=row, col=col)
        fig.add_trace(go.Scatter(
            x=[optimum["search_params"][pair["y_name"]]],
            y=[optimum["search_params"][pair["x_name"]]],
            mode="markers",
            marker={
                "symbol": "star",
                "size": 16,
                "color": "black",
                "line": {"color": "white", "width": 1.5},
            },
            hovertemplate=(
                "GP optimum candidate<br>"
                f"{pair['y_label']}: %{{x:.4g}}<br>"
                f"{pair['x_label']}: %{{y:.4g}}<extra></extra>"
            ),
            showlegend=False,
            meta=trace_meta(pair["id"], "optimum"),
        ), row=row, col=col)
        fig.add_trace(go.Scatter(
            x=[less_informative["search_params"][pair["y_name"]]],
            y=[less_informative["search_params"][pair["x_name"]]],
            mode="markers",
            marker={
                "symbol": "diamond",
                "size": 11,
                "color": "#00B8D9",
                "line": {"color": "white", "width": 1.5},
            },
            hovertemplate=(
                "Less-informative candidate<br>"
                f"{pair['y_label']}: %{{x:.4g}}<br>"
                f"{pair['x_label']}: %{{y:.4g}}<extra></extra>"
            ),
            showlegend=False,
            meta=trace_meta(pair["id"], "less_informative"),
        ), row=row, col=col)
        fig.update_xaxes(
            visible=True,
            type=specs[n_params - col].scale,
            showticklabels=row + col == n_params,
            ticks="outside" if row + col == n_params else "",
            showgrid=False,
            zeroline=False,
            showline=True,
            mirror=True,
            row=row,
            col=col,
        )
        fig.update_yaxes(
            visible=True,
            type=specs[row - 1].scale,
            showticklabels=col == 1,
            ticks="outside" if col == 1 else "",
            showgrid=False,
            zeroline=False,
            showline=True,
            mirror=True,
            row=row,
            col=col,
        )

    for index in range(grid_size):
        col_subplot = fig.get_subplot(1, index + 1)
        row_subplot = fig.get_subplot(index + 1, 1)
        fig.add_annotation(
            x=sum(col_subplot.xaxis.domain) / 2,
            y=1.012,
            xref="paper",
            yref="paper",
            text=specs[n_params - index - 1].label,
            textangle=0,
            showarrow=False,
            xanchor="center",
            yanchor="bottom",
            font={"size": 11},
        )
        fig.add_annotation(
            x=-0.018,
            y=sum(row_subplot.yaxis.domain) / 2,
            xref="paper",
            yref="paper",
            text=specs[index].label,
            showarrow=False,
            xanchor="right",
            font={"size": 12},
        )

    fig.add_annotation(
        x=0.5,
        y=1.095,
        xref="paper",
        yref="paper",
        text="★ GP optimum &nbsp;&nbsp; ◆ less-informative candidate",
        showarrow=False,
        font={"size": 13},
    )
    fig.add_annotation(
        x=0.5,
        y=1.14,
        xref="paper",
        yref="paper",
        text="Pairwise GP posterior trends",
        showarrow=False,
        font={"size": 20},
    )
    fig.update_layout(
        width=figure_size,
        height=figure_size,
        margin={"t": 150, "b": 100, "l": 150, "r": 100},
        showlegend=False,
        coloraxis={
            "colorscale": "Viridis",
            "cmin": report["zmin"],
            "cmax": report["zmax"],
            "colorbar": {
                "title": "predicted accuracy",
                "x": 1.02,
                "xanchor": "left",
                "y": 0.5,
                "len": 0.55,
            },
        },
    )
    return fig


def make_accuracy_figure(report: dict, width: int) -> go.Figure:
    mode = report["mode"]
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("Best observed accuracy", "Observed accuracy and estimated noise"),
        horizontal_spacing=0.13,
    )

    core = [record for record in report["records"] if record["source"] in ("sobol", "bo")]
    running_best = -math.inf
    best = []
    for record in core:
        running_best = max(running_best, record["objective_mean"])
        best.append(running_best)
    fig.add_trace(go.Scatter(
        x=list(range(1, len(core) + 1)), y=best, mode="lines+markers",
        line={"color": MODE_COLORS[mode]}, showlegend=False,
        hovertemplate="evaluation %{x}<br>best accuracy=%{y:.5f}<extra></extra>",
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=list(range(1, len(report["records"]) + 1)),
        y=[record["objective_mean"] for record in report["records"]],
        error_y={"type": "data", "array": [
            math.sqrt(record["objective_variance"]) for record in report["records"]
        ], "visible": True},
        text=[record["source"] for record in report["records"]], mode="markers",
        marker={"color": MODE_COLORS[mode], "size": 7}, showlegend=False,
        hovertemplate="evaluation %{x}<br>accuracy=%{y:.5f}<br>%{text}<extra></extra>",
    ), row=1, col=2)
    fig.update_xaxes(title_text="BO evaluation", row=1, col=1)
    fig.update_xaxes(title_text="evaluation", row=1, col=2)
    fig.update_yaxes(title_text="accuracy", row=1, col=1)
    fig.update_yaxes(title_text="accuracy", row=1, col=2)
    fig.update_layout(
        width=width,
        height=460,
        margin={"t": 75, "b": 70, "l": 80, "r": 45},
        showlegend=False,
    )
    return fig


def make_dataset_comparison_figure(experiment: dict) -> go.Figure:
    datasets = experiment["datasets"]
    modes = list(datasets[0]["modes"])
    fig = go.Figure()
    for mode in modes:
        means = []
        errors = []
        hover = []
        for dataset in datasets:
            confirmation = dataset["modes"][mode]["summary"][
                "predicted_optimum"
            ]["confirmation"]
            means.append(confirmation["objective_mean"])
            errors.append(math.sqrt(confirmation["objective_variance"]))
            hover.append(dataset["seq_file"])
        fig.add_trace(go.Bar(
            name=mode,
            x=[dataset["label"] for dataset in datasets],
            y=means,
            error_y={"type": "data", "array": errors, "visible": True},
            marker_color=MODE_COLORS[mode],
            customdata=hover,
            hovertemplate=(
                "%{customdata}<br>mode=" + mode
                + "<br>actual optimum mean=%{y:.5f}<extra></extra>"
            ),
        ))
    fig.update_layout(
        title="Actual optimum imputation accuracy by sequence file",
        barmode="group",
        width=max(900, min(1500, 150 * len(datasets) + 350)),
        height=520,
        margin={"t": 90, "b": 150, "l": 85, "r": 40},
        xaxis={"title": "sequence file", "tickangle": -25},
        yaxis={"title": "DFCP imputation accuracy"},
        legend={"title": "DFCP mode", "orientation": "h", "y": 1.08},
    )
    return fig


def get_sidebar_context(reports: list[dict]) -> dict:
    context = {"datasets": [], "modes": [], "reports": {}}
    for report in reports:
        dataset = next(
            (
                candidate for candidate in context["datasets"]
                if candidate["id"] == report["dataset_id"]
            ),
            None,
        )
        if dataset is None:
            dataset = {
                "id": report["dataset_id"],
                "label": report["dataset_label"],
                "file": report["seq_file"],
            }
            context["datasets"].append(dataset)
        if report["mode"] not in context["modes"]:
            context["modes"].append(report["mode"])
        summary = report["summary"]
        context["reports"][f'{report["dataset_id"]}|{report["mode"]}'] = {
            "pairs": [{
                "id": pair["id"],
                "label": f"{pair['x_label']} × {pair['y_label']}",
                "x_label": pair["x_label"],
                "y_label": pair["y_label"],
                "x_scale": next(
                    spec.scale for spec in report["specs"]
                    if spec.name == pair["x_name"]
                ),
                "y_scale": next(
                    spec.scale for spec in report["specs"]
                    if spec.name == pair["y_name"]
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
const context = __REPORT_CONTEXT__;
const report = context.reports;
let selectedDataset = context.datasets[0].id;
let selectedMode = context.modes[0];
let selectedPair = report[`${selectedDataset}|${selectedMode}`].pairs[0].id;
let reportWidth;

function selectedReport() {
    return report[`${selectedDataset}|${selectedMode}`];
}

function selectedSuffix() {
    return `${selectedDataset}-${selectedMode}`;
}

function updatePairOptions() {
    const select = document.getElementById('pair-select');
    select.innerHTML = selectedReport().pairs.map(pair =>
        `<option value="${pair.id}">${pair.label}</option>`
    ).join('');
    selectedPair = selectedReport().pairs[0].id;
    select.value = selectedPair;
}

function updateOverview() {
    const plot = document.getElementById(`tune-overview-${selectedSuffix()}`);
    if (!plot) return;
    const pair = selectedReport().pairs.find(candidate => candidate.id === selectedPair);
    const grid = document.getElementById(`tune-grid-${selectedSuffix()}`);
    const transpose = values => {
        if (!values || values.length === 0) return values;
        return values[0].map((_, index) => values.map(row => row[index]));
    };
    const traces = grid.data.filter(trace =>
        trace.meta?.kind === 'grid_surface' && trace.meta.pair === selectedPair
    ).map(source => {
        const trace = {
            ...source,
            x: source.y,
            y: source.x,
            xaxis: 'x2',
            yaxis: 'y2',
            meta: {kind: 'overview_surface', role: source.meta.role},
        };
        if (source.z) trace.z = transpose(source.z);
        if (source.customdata) trace.customdata = transpose(source.customdata);
        if (source.meta.role === 'heatmap') {
            trace.hovertemplate = `${pair.x_label}: %{x:.4g}<br>` +
                `${pair.y_label}: %{y:.4g}<br>predicted accuracy: %{z:.5f}<br>` +
                'posterior SD: %{customdata:.5f}<extra></extra>';
        } else if (source.meta.role === 'contour') {
            trace.contours = {...source.contours, showlabels: true};
            trace.textfont = {...source.textfont, size: 11};
        } else if (source.meta.role === 'observed') {
            trace.marker = {...source.marker, size: 5, opacity: 1.0};
        } else if (source.meta.role === 'optimum') {
            trace.hovertemplate = 'GP optimum candidate<extra></extra>';
        } else if (source.meta.role === 'less_informative') {
            trace.hovertemplate = 'Less-informative candidate<extra></extra>';
        }
        return trace;
    });
    Plotly.react(plot, [plot.data[0], ...traces], plot.layout, plot._context);
    Plotly.relayout(plot, {
        'xaxis2.title.text': pair.x_label,
        'xaxis2.type': pair.x_scale,
        'xaxis2.autorange': true,
        'yaxis2.title.text': pair.y_label,
        'yaxis2.type': pair.y_scale,
        'yaxis2.autorange': true,
    });
}

function updateSummary() {
    const mode = selectedReport();
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

function updatePlot() {
    document.querySelectorAll('.dataset-mode-report').forEach(section => {
        section.hidden = (
            section.dataset.dataset !== selectedDataset ||
            section.dataset.mode !== selectedMode
        );
    });
    updateSummary();
    requestAnimationFrame(() => {
        resizePlots();
        updateOverview();
    });
}

function resizePlots() {
    const section = document.querySelector(
        `.dataset-mode-report[data-dataset="${selectedDataset}"][data-mode="${selectedMode}"]`
    );
    if (!section || section.hidden) return;
    const plots = {
        overview: document.getElementById(`tune-overview-${selectedSuffix()}`),
        grid: document.getElementById(`tune-grid-${selectedSuffix()}`),
        accuracy: document.getElementById(`tune-accuracy-${selectedSuffix()}`),
    };
    const availableWidth = Math.max(760, Math.min(section.clientWidth * 0.94, 1360));
    reportWidth = reportWidth === undefined
        ? availableWidth
        : Math.min(reportWidth, availableWidth);
    const targetWidth = reportWidth;
    plots.overview.parentElement.style.width = `${targetWidth}px`;
    plots.grid.parentElement.style.width = `${targetWidth}px`;
    plots.grid.parentElement.style.height = `${targetWidth}px`;
    plots.accuracy.parentElement.style.width = `${targetWidth}px`;
    Plotly.relayout(plots.overview, {width: targetWidth});
    Plotly.relayout(plots.grid, {width: targetWidth, height: targetWidth});
    Plotly.relayout(plots.accuracy, {width: targetWidth});
    const comparison = document.getElementById('tune-dataset-comparison');
    if (comparison) {
        comparison.parentElement.style.width = `${targetWidth}px`;
        Plotly.relayout(comparison, {width: targetWidth});
    }
}

const sidebar = document.createElement('aside');
sidebar.id = 'tune-sidebar';
sidebar.innerHTML = `
    <h2>Hyperparameter report</h2>
    <section><h3>Haplotype dataset</h3>
        <select id="dataset-select">${context.datasets.map(dataset =>
            `<option value="${dataset.id}">${dataset.label}</option>`
        ).join('')}</select>
        <div class="tune-dataset-file" id="dataset-file">${context.datasets[0].file}</div>
    </section>
    <section><h3>DFCP mode</h3><div class="tune-buttons" id="mode-controls">
        ${context.modes.map((mode, index) =>
            `<button class="${index === 0 ? 'active' : ''}" data-mode="${mode}">${mode}</button>`
        ).join('')}
    </div></section>
    <section><h3>Large contour pair</h3><select id="pair-select"></select></section>
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
        <p>Each heatmap is the GP posterior mean with unshown parameters fixed at the
        less-informative candidate. White lines show equal predicted accuracy.</p>
        <p>The black star is the GP optimum, the cyan diamond is the less-informative
        candidate, and dots are evaluated configurations projected into each cell.</p>
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
    #pair-select, #dataset-select { box-sizing: border-box; width: 100%; padding: 7px 5px;
        border: 1px solid #aeb5c0; border-radius: 4px; background: white;
        color: #2a3f5f; }
    .tune-dataset-file { margin-top: 6px; overflow-wrap: anywhere;
        color: #697386; font-size: 11px; line-height: 1.3; }
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
    .dataset-mode-report, .comparison-report { width: calc(100% - 350px);
        margin-left: 350px; overflow-x: auto; }
    .dataset-mode-report[hidden] { display: none; }
    .plot-stack { width: 100%; margin: 0 auto; }
    .plot-section { margin-bottom: 44px; }
    .plot-section-grid { margin-bottom: 44px; }
    @media (max-width: 850px) { #tune-sidebar { width: 280px; padding: 12px; }
        .dataset-mode-report, .comparison-report {
            width: calc(100% - 280px); margin-left: 280px;
        } }`;
document.head.appendChild(style);

const datasetSelect = document.getElementById('dataset-select');
function applySelectedDataset() {
    selectedDataset = datasetSelect.value;
    document.getElementById('dataset-file').textContent = context.datasets.find(
        dataset => dataset.id === selectedDataset
    ).file;
    updatePairOptions();
    updatePlot();
}
datasetSelect.addEventListener('input', applySelectedDataset);
datasetSelect.addEventListener('change', applySelectedDataset);

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
const pairSelect = document.getElementById('pair-select');
function applySelectedPair() {
    selectedPair = pairSelect.value;
    updateOverview();
}
pairSelect.addEventListener('input', applySelectedPair);
pairSelect.addEventListener('change', applySelectedPair);
updatePairOptions();
updatePlot();
window.addEventListener('resize', () => requestAnimationFrame(resizePlots));
"""
    return template.replace("__REPORT_CONTEXT__", json.dumps(context))


def write_report(experiment: dict, output: Path) -> None:
    plotly_asset = ensure_plotly_asset(output)
    reports = [
        get_mode_report(dataset_index, dataset, mode, mode_data)
        for dataset_index, dataset in enumerate(experiment["datasets"])
        for mode, mode_data in dataset["modes"].items()
    ]
    fragments = []
    for index, report in enumerate(reports):
        width = PAIR_CELL_SIZE * len(report["specs"]) + 250
        figures = (
            ("overview", make_overview_figure(report, width)),
            ("grid", make_pair_grid_figure(report)),
            ("accuracy", make_accuracy_figure(report, width)),
        )
        plots = []
        suffix = f'{report["dataset_id"]}-{report["mode"]}'
        for name, figure in figures:
            plot = pio.to_html(
                figure,
                full_html=False,
                include_plotlyjs=False,
                config={"responsive": False},
                div_id=f"tune-{name}-{suffix}",
            )
            plots.append(f'<div class="plot-section plot-section-{name}">{plot}</div>')
        hidden = "" if index == 0 else " hidden"
        fragments.append(
            f'<section class="dataset-mode-report" '
            f'data-dataset="{report["dataset_id"]}" '
            f'data-mode="{report["mode"]}"{hidden}>'
            f'<div class="plot-stack">{"".join(plots)}</div></section>'
        )

    comparison = ""
    if len(experiment["datasets"]) > 1:
        comparison_plot = pio.to_html(
            make_dataset_comparison_figure(experiment),
            full_html=False,
            include_plotlyjs=False,
            config={"responsive": False},
            div_id="tune-dataset-comparison",
        )
        comparison = (
            '<section class="comparison-report"><div class="plot-stack">'
            f'<div class="plot-section plot-section-comparison">{comparison_plot}</div>'
            "</div></section>"
        )

    page = (
        '<!doctype html>\n<html><head><meta charset="utf-8">'
        "<title>DFCP hyperparameter report</title>"
        f'<script src="{plotly_asset}"></script></head><body>\n'
        + "\n".join(fragments)
        + comparison
        + "\n<script>\n"
        + get_page_script(get_sidebar_context(reports))
        + "\n</script>\n</body></html>\n"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page)


def save_experiment(experiment: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(experiment, indent=2, allow_nan=False) + "\n")


def load_experiment(path: Path) -> dict:
    experiment = json.loads(path.read_text())
    if not isinstance(experiment, dict):
        raise ValueError(f"{path} is not a DFCP tuning result")
    if isinstance(experiment.get("datasets"), list):
        if not experiment["datasets"]:
            raise ValueError(f"{path} contains no tuned datasets")
        for dataset in experiment["datasets"]:
            if not isinstance(dataset.get("modes"), dict) or not dataset["modes"]:
                raise ValueError(f"{path} contains a dataset with no tuned modes")
        return experiment

    # Version 2 averaged all configured files into one objective. A single-file
    # result is losslessly representable in the per-dataset version 3 schema.
    if isinstance(experiment.get("modes"), dict) and experiment["modes"]:
        seq_files = experiment.get("config", {}).get("seq_files", [])
        if len(seq_files) != 1:
            raise ValueError(
                f"{path} averages {len(seq_files)} sequence files and cannot be "
                "visualized as independent tuning results; rerun tuning"
            )
        experiment = {
            **experiment,
            "version": 3,
            "datasets": [{
                "seq_file": seq_files[0],
                "label": get_seq_label(Path(seq_files[0])),
                "modes": experiment["modes"],
            }],
        }
        experiment.pop("modes")
        return experiment
    raise ValueError(f"{path} is not a DFCP tuning result")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Tune DFCP priors independently per sequence file with BoTorch, save "
            "JSON results, and write a Plotly report."
        )
    )
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--seq_files", type=Path, nargs="+")
    inputs.add_argument(
        "--seq_dir", type=Path,
        help=(
            "independently tune every haps*.txt* file in this directory; by "
            f"default, use only the clean file {DEFAULT_SEQ_FILE}"
        ),
    )
    inputs.add_argument(
        "--load", type=Path, metavar="JSON",
        help="skip tuning and build the Plotly report from saved JSON results",
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
    parser.add_argument("--json", type=Path, default=Path("output/tune.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/tune.html"))
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.load is not None:
        if not args.load.is_file():
            raise ValueError(f"tuning result does not exist: {args.load}")
        return
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
    if args.load is not None:
        experiment = load_experiment(args.load)
        torch.manual_seed(experiment.get("config", {}).get("seed", 0))
        write_report(experiment, args.output)
        print(f"loaded {args.load}")
        print(f"wrote {args.output}")
        return

    torch.manual_seed(args.seed)
    seq_files = resolve_seq_files(args)
    build_dfcp()
    experiment = {
        "version": 3,
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
            "dataset_handling": "independent optimization per sequence file",
            "parameterization": (
                "Gamma mean/shape, discount Beta v1/v2, mismatch Beta mean/concentration"
            ),
        },
        "datasets": [],
    }
    for dataset_index, seq_file in enumerate(seq_files):
        dataset = {
            "seq_file": str(seq_file),
            "label": get_seq_label(seq_file),
            "modes": {},
        }
        for mode in args.modes:
            dataset["modes"][mode] = tune_mode(
                mode, seq_file, dataset_index, args
            )
        experiment["datasets"].append(dataset)

    save_experiment(experiment, args.json)
    write_report(load_experiment(args.json), args.output)
    print(f"wrote {args.json}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
