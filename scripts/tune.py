"""Independent noisy Bayesian optimization of DFCP priors.

Each sequence file and DFCP mode gets its own optimization. A run evaluates the
current defaults, a scrambled Sobol design, sequential qLogNEI suggestions, and
then replicated confirmations of the GP optimum and a near-optimal broad-prior
candidate. Results are JSON only; scripts/tune_viz.py owns all visualization.
"""

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import torch
from botorch.acquisition.logei import qLogNoisyExpectedImprovement
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms.outcome import Standardize
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from gpytorch.mlls import ExactMarginalLogLikelihood

from dfcp import build_dfcp, run_dfcp


DEFAULT_SEQ_DIR = Path("data/examples/simulated/SIM1_LEN500_NHAPS100")
DEFAULT_SEQ_FILE = DEFAULT_SEQ_DIR / "haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt"
PBWT_MATCH_LEN = 20
NEAR_BEST_TOLERANCE = 0.0025
RECOMMENDATION_POINTS = 8192
MC_SAMPLES = 128
ACQUISITION_RESTARTS = 10
ACQUISITION_RAW_SAMPLES = 256


@dataclass(frozen=True)
class ParameterSpec:
    name: str
    lower: float
    upper: float
    scale: str
    group: str

    def from_unit(self, value: float) -> float:
        if self.scale == "log":
            return math.exp(
                math.log(self.lower) + value * (math.log(self.upper) - math.log(self.lower))
            )
        return self.lower + value * (self.upper - self.lower)

    def to_unit(self, value: float) -> float:
        if self.scale == "log":
            return (math.log(value) - math.log(self.lower)) / (
                math.log(self.upper) - math.log(self.lower)
            )
        return (value - self.lower) / (self.upper - self.lower)


COMMON_PARAMETERS = (
    ParameterSpec("alpha_shape", 0.5, 10.0, "log", "alpha"),
    ParameterSpec("alpha_mean", 0.001, 10.0, "log", "alpha"),

    ParameterSpec("d_mean", 0.1, 0.9, "linear", "d"),
    ParameterSpec("d_conc", 1.0, 20.0, "log", "d"),

    ParameterSpec("gamma_shape", 1.0, 10.0, "log", "gamma"),
    ParameterSpec("gamma_mean", 0.001, 10.0, "log", "gamma"),
)
NOISY_PARAMETERS = (
    ParameterSpec("eps_mean", 0.001, 0.2, "log", "eps"),
    ParameterSpec("eps_conc", 1.0, 10.0, "log", "eps"),
)


def parameter_specs(mode: str) -> tuple[ParameterSpec, ...]:
    return COMMON_PARAMETERS + (NOISY_PARAMETERS if mode == "noisy" else ())


def search_params_from_unit(
    point: list[float], specs: tuple[ParameterSpec, ...]
) -> dict[str, float]:
    return {spec.name: spec.from_unit(value) for spec, value in zip(specs, point)}


def cpp_params(search_params: dict[str, float], mode: str) -> dict[str, float]:
    """Convert Gamma shape/mean and Beta mean/concentration to C++ parameters."""
    d_mean = search_params["d_mean"]
    d_conc = search_params["d_conc"]
    result = {
        "tau_1": search_params["alpha_shape"],
        "tau_2": search_params["alpha_shape"] / search_params["alpha_mean"],
        "v_1": d_mean * d_conc,
        "v_2": (1.0 - d_mean) * d_conc,
        "phi_1": search_params["gamma_shape"],
        "phi_2": search_params["gamma_shape"] / search_params["gamma_mean"],
    }
    if mode == "noisy":
        eps_mean = search_params["eps_mean"]
        eps_conc = search_params["eps_conc"]
        result["lambda_1"] = eps_mean * eps_conc
        result["lambda_2"] = (1.0 - eps_mean) * eps_conc
    return result


def default_unit_point(specs: tuple[ParameterSpec, ...]) -> list[float]:
    defaults = {
        "alpha_shape": 1.0,
        "alpha_mean": 1.0,
        "d_mean": 0.5,
        "d_conc": 2.0,
        "gamma_shape": 2.0,
        "gamma_mean": 1.0,
        "eps_mean": 0.1 / 1.1,
        "eps_conc": 1.1,
    }
    return [spec.to_unit(defaults[spec.name]) for spec in specs]


def evaluate_point(
    point: list[float],
    specs: tuple[ParameterSpec, ...],
    mode: str,
    seq_file: Path,
    mask: float,
    val: float,
    replicates: int,
    source: str,
    iteration: int,
) -> dict:
    """Estimate one candidate's mean accuracy and the noise variance of that mean.

    Every replicate reruns DFCP on the same sequence file. DFCP independently
    samples validation sequences and masked alleles, so replicate accuracy varies.
    The GP observation variance is the larger of the empirical variance of the
    replicate mean and the average binomial sampling variance of that mean.
    """
    search_params = search_params_from_unit(point, specs)
    dfcp_params = cpp_params(search_params, mode)
    runs = []
    for replicate in range(replicates):
        output = run_dfcp(
            seq_file,
            retries=1,
            pbwt_init=1,
            pbwt_match_len=PBWT_MATCH_LEN,
            mask=mask,
            val=val,
            soft=int(mode == "soft"),
            noisy=int(mode == "noisy"),
            **dfcp_params,
        )
        runs.append(
            {
                "replicate": replicate,
                "accuracy": float(output["dfcp_impute_acc"]),
                "n_masked_alleles": int(output["n_masked_alleles"]),
            }
        )

    scores = [run["accuracy"] for run in runs]
    objective_mean = sum(scores) / replicates
    sample_variance = sum((score - objective_mean) ** 2 for score in scores) / (replicates - 1)
    variance_of_mean = sample_variance / replicates
    binomial_variance_of_mean = sum(
        run["accuracy"] * (1.0 - run["accuracy"]) / run["n_masked_alleles"] for run in runs
    ) / (replicates**2)
    objective_variance = max(variance_of_mean, binomial_variance_of_mean, 1e-8)

    record = {
        "x": point,
        "search_params": search_params,
        "cpp_params": dfcp_params,
        "objective_mean": objective_mean,
        "objective_variance": objective_variance,
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
    train_y = torch.tensor([[record["objective_mean"]] for record in records], dtype=torch.double)
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


def propose_next_point(records: list[dict]) -> list[float]:
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
        num_restarts=ACQUISITION_RESTARTS,
        raw_samples=ACQUISITION_RAW_SAMPLES,
    )
    return candidate.squeeze(0).tolist()


def prior_information_penalty(point: torch.Tensor, specs: tuple[ParameterSpec, ...]) -> float:
    """Rank near-optimal candidates by distance from broad reference priors."""
    params = search_params_from_unit(point.tolist(), specs)

    def beta_penalty(mean: float, conc: float) -> float:
        # Beta(1,1), with mean 0.5 and concentration 2, is the broad reference.
        return abs(math.log(mean / (1.0 - mean))) + abs(math.log(conc / 2.0))

    penalty = math.log(params["alpha_shape"]) + math.log(params["gamma_shape"])
    penalty += beta_penalty(params["d_mean"], params["d_conc"])
    if "eps_mean" in params:
        penalty += beta_penalty(params["eps_mean"], params["eps_conc"])
    return penalty


def select_candidates(
    model: SingleTaskGP,
    specs: tuple[ParameterSpec, ...],
    mode: str,
    seed: int,
) -> tuple[dict, dict]:
    candidates = (
        torch.quasirandom.SobolEngine(len(specs), scramble=True, seed=seed)
        .draw(RECOMMENDATION_POINTS)
        .double()
    )
    with torch.no_grad():
        posterior_means = model.posterior(candidates).mean.squeeze(-1)

    optimum_index = int(torch.argmax(posterior_means))
    best_prediction = float(posterior_means[optimum_index])
    near_best = torch.nonzero(posterior_means >= best_prediction - NEAR_BEST_TOLERANCE).squeeze(-1)
    less_informative_index = min(
        near_best.tolist(),
        key=lambda index: prior_information_penalty(candidates[index], specs),
    )

    def describe(index: int) -> dict:
        point = candidates[index].tolist()
        search_params = search_params_from_unit(point, specs)
        return {
            "x": point,
            "posterior_mean": float(posterior_means[index]),
            "search_params": search_params,
            "cpp_params": cpp_params(search_params, mode),
        }

    return describe(optimum_index), describe(less_informative_index)


def ard_summary(
    model: SingleTaskGP, specs: tuple[ParameterSpec, ...]
) -> tuple[dict[str, float], dict[str, float]]:
    values = model.covar_module.lengthscale.detach().reshape(-1).tolist()
    lengthscales = {spec.name: value for spec, value in zip(specs, values)}
    inverse = {name: 1.0 / value for name, value in lengthscales.items()}
    total = sum(inverse.values())
    normalized_inverse = {name: value / total for name, value in inverse.items()}
    return lengthscales, normalized_inverse


def tune_mode(
    mode: str,
    seq_file: Path,
    dataset_index: int,
    args: argparse.Namespace,
) -> dict:
    """Run Sobol initialization, noisy BO, and candidate confirmation for one mode."""
    specs = parameter_specs(mode)
    mode_index = args.modes.index(mode)
    seed = args.seed + 100_000 * dataset_index + 1_000 * mode_index

    current_default = evaluate_point(
        default_unit_point(specs),
        specs,
        mode,
        seq_file,
        args.mask,
        args.val,
        args.confirmation_replicates,
        "current_default",
        -1,
    )

    records = [current_default]
    sobol = torch.quasirandom.SobolEngine(len(specs), scramble=True, seed=seed)
    for iteration, point in enumerate(sobol.draw(args.initial_points).double()):
        records.append(
            evaluate_point(
                point.tolist(),
                specs,
                mode,
                seq_file,
                args.mask,
                args.val,
                args.replicates,
                "sobol",
                iteration,
            )
        )

    for iteration in range(args.iterations):
        point = propose_next_point(records)
        records.append(
            evaluate_point(
                point,
                specs,
                mode,
                seq_file,
                args.mask,
                args.val,
                args.replicates,
                "bo",
                iteration,
            )
        )

    selection_model = fit_surrogate(records)
    predicted_optimum, less_informative = select_candidates(
        selection_model, specs, mode, seed + 10_000
    )
    search_records = list(records)
    predicted_optimum["confirmation"] = evaluate_point(
        predicted_optimum["x"],
        specs,
        mode,
        seq_file,
        args.mask,
        args.val,
        args.confirmation_replicates,
        "confirm_predicted_optimum",
        len(records),
    )
    less_informative["confirmation"] = evaluate_point(
        less_informative["x"],
        specs,
        mode,
        seq_file,
        args.mask,
        args.val,
        args.confirmation_replicates,
        "confirm_less_informative",
        len(records) + 1,
    )
    records.extend((predicted_optimum["confirmation"], less_informative["confirmation"]))

    final_model = fit_surrogate(records)
    lengthscales, normalized_inverse = ard_summary(final_model, specs)
    return {
        "parameter_specs": [asdict(spec) for spec in specs],
        "records": records,
        "summary": {
            "current_default": current_default,
            "best_observed": max(search_records, key=lambda record: record["objective_mean"]),
            "predicted_optimum": predicted_optimum,
            "less_informative": less_informative,
            "ard_lengthscales": lengthscales,
            "normalized_inverse_lengthscales": normalized_inverse,
        },
    }


def resolve_seq_files(args: argparse.Namespace) -> list[Path]:
    if args.seq_files:
        files = args.seq_files
    elif args.seq_dir:
        files = sorted(args.seq_dir.glob("haps*.txt*"))
    else:
        files = [DEFAULT_SEQ_FILE]
    if not files:
        raise ValueError("no haplotype files were found")
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(f"missing sequence file: {path}")
    return files


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tune DFCP priors independently for each sequence file."
    )
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--seq_files", type=Path, nargs="+")
    inputs.add_argument(
        "--seq_dir",
        type=Path,
        help="independently tune every haps*.txt* file in this directory",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=("hard", "noisy", "soft"),
        default=["hard", "noisy", "soft"],
    )
    parser.add_argument("--mask", type=float, default=0.2)
    parser.add_argument("--val", type=float, default=0.2)
    parser.add_argument("--initial_points", type=int, default=16)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--confirmation_replicates", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=Path("output/tune.json"))
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
    if args.replicates < 2 or args.confirmation_replicates < 2:
        raise ValueError("replicate counts must be at least 2")


def main() -> None:
    args = parse_args()
    validate_args(args)
    torch.manual_seed(args.seed)
    seq_files = resolve_seq_files(args)
    build_dfcp()

    experiment = {
        "config": {
            "mask": args.mask,
            "val": args.val,
            "pbwt_match_len": PBWT_MATCH_LEN,
            "initial_points": args.initial_points,
            "iterations": args.iterations,
            "replicates": args.replicates,
            "confirmation_replicates": args.confirmation_replicates,
            "seed": args.seed,
        },
        "datasets": [],
    }
    for dataset_index, seq_file in enumerate(seq_files):
        dataset = {"seq_file": str(seq_file), "modes": {}}
        for mode in args.modes:
            dataset["modes"][mode] = tune_mode(mode, seq_file, dataset_index, args)
        experiment["datasets"].append(dataset)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(experiment, indent=2, allow_nan=False) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
