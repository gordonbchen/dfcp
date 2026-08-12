"""Tune DFCP priors independently for each sequence file and mode.

The tuner evaluates the current defaults once, builds broad coverage with a
scrambled Sobol design, and then performs sequential Bayesian optimization with
qLogNoisyExpectedImprovement. Every hyperparameter location gets one DFCP run.
A homoskedastic observation-noise variance is inferred with the other GP
hyperparameters instead of being estimated from fixed replications.

After the main search, adaptive recommendation rounds repeatedly optimize the
GP posterior mean and a constrained less-informative objective, evaluate both
recommendations once, and refit. The final recommendations are recomputed from
all observations and selected from locations DFCP actually evaluated. This
spends repeated evaluations only where the updated GP continues to recommend
them. Results are JSON only; scripts/tune_viz.py owns the interactive report.
"""

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import torch
from botorch.acquisition import PosteriorMean
from botorch.acquisition.acquisition import AcquisitionFunction
from botorch.acquisition.logei import qLogNoisyExpectedImprovement
from botorch.exceptions.errors import CandidateGenerationError
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms.outcome import Standardize
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.utils.transforms import t_batch_mode_transform
from gpytorch.mlls import ExactMarginalLogLikelihood

from dfcp import build_dfcp, run_dfcp
from seq_file_name import get_seq_sort_key


DEFAULT_SEQ_DIR = Path("data/examples/simulated/SIM1_LEN500_NHAPS100")
DEFAULT_SEQ_FILE = DEFAULT_SEQ_DIR / "haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt"
PBWT_MATCH_LEN = 20
NEAR_BEST_TOLERANCE = 0.0025

# These control inexpensive numerical optimization of GP-based objectives, not
# the number of DFCP evaluations.
ACQUISITION_MC_SAMPLES = 256
ACQUISITION_RESTARTS = 20
ACQUISITION_RAW_SAMPLES = 1024

RECOMMENDATION_RESTARTS = 20
RECOMMENDATION_RAW_SAMPLES = 2048
CONSTRAINED_RAW_SAMPLES = 16384
CONSTRAINT_START_MARGIN = 1e-5


@dataclass(frozen=True)
class ParameterSpec:
    name: str
    lower: float
    upper: float
    scale: str
    group: str

    def from_unit(self, value: float) -> float:
        if self.scale == "log":
            log_lower = math.log(self.lower)
            return math.exp(log_lower + value * (math.log(self.upper) - log_lower))
        return self.lower + value * (self.upper - self.lower)

    def to_unit(self, value: float) -> float:
        if self.scale == "log":
            return (math.log(value) - math.log(self.lower)) / (
                math.log(self.upper) - math.log(self.lower)
            )
        return (value - self.lower) / (self.upper - self.lower)

    def from_unit_tensor(self, value: torch.Tensor) -> torch.Tensor:
        if self.scale == "log":
            log_lower = math.log(self.lower)
            return torch.exp(log_lower + value * (math.log(self.upper) - log_lower))
        return self.lower + value * (self.upper - self.lower)


COMMON_PARAMETERS = (
    ParameterSpec("alpha_shape", 0.1, 10.0, "log", "alpha"),
    ParameterSpec("alpha_mean", 0.001, 10.0, "log", "alpha"),

    ParameterSpec("d_mean", 0.1, 0.9, "linear", "d"),
    ParameterSpec("d_conc", 1.0, 20.0, "log", "d"),

    ParameterSpec("gamma_shape", 0.1, 10.0, "log", "gamma"),
    ParameterSpec("gamma_mean", 0.0001, 5.0, "log", "gamma"),
)
NOISY_PARAMETERS = (
    ParameterSpec("eps_mean", 0.001, 0.2, "log", "eps"),
    ParameterSpec("eps_conc", 1.0, 10.0, "log", "eps"),
)


def parameter_specs(mode: str) -> tuple[ParameterSpec, ...]:
    return COMMON_PARAMETERS + (NOISY_PARAMETERS if mode == "noisy" else ())


def unit_bounds(dimension: int) -> torch.Tensor:
    lower = torch.zeros(dimension, dtype=torch.double)
    upper = torch.ones(dimension, dtype=torch.double)
    return torch.stack((lower, upper))


def search_params_from_unit(
    point: list[float], specs: tuple[ParameterSpec, ...]
) -> dict[str, float]:
    return {spec.name: spec.from_unit(value) for spec, value in zip(specs, point)}


def tensor_params_from_unit(
    points: torch.Tensor, specs: tuple[ParameterSpec, ...]
) -> dict[str, torch.Tensor]:
    return {
        spec.name: spec.from_unit_tensor(points[..., index])
        for index, spec in enumerate(specs)
    }


def cpp_params(search_params: dict[str, float], mode: str) -> dict[str, float]:
    """Convert Gamma shape/mean and Beta mean/concentration to C++ parameters."""
    alpha_shape = search_params["alpha_shape"]
    alpha_mean = search_params["alpha_mean"]
    d_mean = search_params["d_mean"]
    d_conc = search_params["d_conc"]
    gamma_shape = search_params["gamma_shape"]
    gamma_mean = search_params["gamma_mean"]
    result = {
        "tau_1": alpha_shape, "tau_2": alpha_shape / alpha_mean,
        "v_1": d_mean * d_conc, "v_2": (1.0 - d_mean) * d_conc,
        "phi_1": gamma_shape, "phi_2": gamma_shape / gamma_mean,
    }
    if mode == "noisy":
        eps_mean = search_params["eps_mean"]
        eps_conc = search_params["eps_conc"]
        result["lambda_1"] = eps_mean * eps_conc
        result["lambda_2"] = (1.0 - eps_mean) * eps_conc
    return result


def default_unit_point(specs: tuple[ParameterSpec, ...]) -> list[float]:
    defaults = {
        "alpha_shape": 1.0, "alpha_mean": 1.0,
        "d_mean": 0.5, "d_conc": 2.0,
        "gamma_shape": 2.0, "gamma_mean": 1.0,
        "eps_mean": 0.1 / 1.1, "eps_conc": 1.1,
    }
    return [spec.to_unit(defaults[spec.name]) for spec in specs]


def evaluate_point(
    point: list[float], specs: tuple[ParameterSpec, ...], mode: str, seq_file: Path,
    mask: float, val: float, source: str, iteration: int,
) -> dict:
    """Run one independently randomized DFCP evaluation."""
    search_params = search_params_from_unit(point, specs)
    dfcp_params = cpp_params(search_params, mode)
    output = run_dfcp(
        seq_file, retries=1,
        pbwt_init=1, pbwt_match_len=PBWT_MATCH_LEN,
        mask=mask, val=val,
        soft=int(mode == "soft"), noisy=int(mode == "noisy"),
        **dfcp_params,
    )
    record = {
        "x": point,
        "search_params": search_params,
        "cpp_params": dfcp_params,
        "accuracy": float(output["viterbi_impute_acc"]),
        "n_masked_ls": int(output["n_masked_ls"]),
        "source": source,
        "iteration": iteration,
    }
    print(
        f"{seq_file.name} {mode:5s} {source:28s} {iteration:3d} "
        f"accuracy={record['accuracy']:.6f}",
        flush=True,
    )
    return record


def fit_surrogate(records: list[dict]) -> SingleTaskGP:
    """Fit an exact GP with one learned, shared observation-noise variance."""
    train_x = torch.tensor([record["x"] for record in records], dtype=torch.double)
    train_y = torch.tensor([[record["accuracy"]] for record in records], dtype=torch.double)
    model = SingleTaskGP(train_x, train_y, outcome_transform=Standardize(m=1))
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))
    return model


def observation_noise_variance(model: SingleTaskGP) -> float:
    """Return learned observation noise in the original accuracy scale."""
    standardized_noise = model.likelihood.noise.detach().reshape(-1)[0]
    outcome_std = model.outcome_transform.stdvs.detach().reshape(-1)[0]
    return float(standardized_noise * outcome_std.square())


def propose_next_point(records: list[dict]) -> list[float]:
    model = fit_surrogate(records)
    train_x = torch.tensor([record["x"] for record in records], dtype=torch.double)
    acquisition = qLogNoisyExpectedImprovement(
        model=model, X_baseline=train_x,
        sampler=SobolQMCNormalSampler(torch.Size([ACQUISITION_MC_SAMPLES])),
    )
    candidate, _ = optimize_acqf(
        acquisition, bounds=unit_bounds(train_x.shape[1]), q=1,
        num_restarts=ACQUISITION_RESTARTS, raw_samples=ACQUISITION_RAW_SAMPLES,
    )
    return candidate.squeeze(0).tolist()


def prior_information_penalty(
    points: torch.Tensor, specs: tuple[ParameterSpec, ...]
) -> torch.Tensor:
    """Score distance from broad Gamma-shape and uniform-Beta references."""
    params = tensor_params_from_unit(points, specs)

    def beta_penalty(mean: torch.Tensor, conc: torch.Tensor) -> torch.Tensor:
        return torch.abs(torch.logit(mean)) + torch.abs(torch.log(conc / 2.0))

    penalty = torch.log(params["alpha_shape"]) + torch.log(params["gamma_shape"])
    penalty = penalty + beta_penalty(params["d_mean"], params["d_conc"])
    if "eps_mean" in params:
        penalty = penalty + beta_penalty(params["eps_mean"], params["eps_conc"])
    return penalty


class NegativePriorInformation(AcquisitionFunction):
    """A differentiable objective for constrained recommendation optimization."""

    def __init__(self, model: SingleTaskGP, specs: tuple[ParameterSpec, ...]):
        super().__init__(model)
        self.specs = specs

    @t_batch_mode_transform(expected_q=1)
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return -prior_information_penalty(x.squeeze(-2), self.specs)


def describe_candidate(
    model: SingleTaskGP, point: torch.Tensor,
    specs: tuple[ParameterSpec, ...], mode: str,
) -> dict:
    with torch.no_grad():
        posterior = model.posterior(point.unsqueeze(0))
        posterior_mean = float(posterior.mean.squeeze())
        posterior_sd = float(posterior.variance.clamp_min(0.0).sqrt().squeeze())
        penalty = float(prior_information_penalty(point.unsqueeze(0), specs).squeeze())
    point_list = point.tolist()
    search_params = search_params_from_unit(point_list, specs)
    return {
        "x": point_list,
        "posterior_mean": posterior_mean,
        "posterior_sd": posterior_sd,
        "information_penalty": penalty,
        "search_params": search_params,
        "cpp_params": cpp_params(search_params, mode),
    }


def select_evaluated_candidates(
    model: SingleTaskGP, records: list[dict],
    specs: tuple[ParameterSpec, ...], mode: str,
) -> tuple[dict, dict]:
    """Select final recommendations from locations that DFCP evaluated."""
    points = torch.tensor([record["x"] for record in records], dtype=torch.double)
    with torch.no_grad():
        means = model.posterior(points).mean.squeeze(-1)
        optimum_index = int(torch.argmax(means))
        threshold = means[optimum_index] - NEAR_BEST_TOLERANCE
        feasible = torch.nonzero(means >= threshold).squeeze(-1)
        penalties = prior_information_penalty(points[feasible], specs)
        less_informative_index = int(feasible[torch.argmin(penalties)])

    def describe(index: int) -> dict:
        candidate = describe_candidate(model, points[index], specs, mode)
        candidate["evaluation_index"] = index
        return candidate

    return describe(optimum_index), describe(less_informative_index)


def optimize_posterior_mean(model: SingleTaskGP, dimension: int) -> torch.Tensor:
    candidate, _ = optimize_acqf(
        PosteriorMean(model), bounds=unit_bounds(dimension), q=1,
        num_restarts=RECOMMENDATION_RESTARTS, raw_samples=RECOMMENDATION_RAW_SAMPLES,
    )
    return candidate.squeeze(0)


def optimize_less_informative(
    model: SingleTaskGP, specs: tuple[ParameterSpec, ...],
    optimum: torch.Tensor, seed: int,
) -> torch.Tensor:
    """Minimize prior information subject to near-optimal posterior accuracy."""
    with torch.no_grad():
        optimum_mean = model.posterior(optimum.unsqueeze(0)).mean.squeeze()
    threshold = optimum_mean - NEAR_BEST_TOLERANCE

    raw_points = (
        torch.quasirandom.SobolEngine(len(specs), scramble=True, seed=seed)
        .draw(CONSTRAINED_RAW_SAMPLES)
        .double()
    )
    raw_points = torch.cat((optimum.unsqueeze(0), raw_points))
    with torch.no_grad():
        raw_means = torch.cat(
            [model.posterior(chunk).mean.squeeze(-1) for chunk in raw_points.split(4096)]
        )
        feasible = raw_points[raw_means >= threshold]
        penalties = prior_information_penalty(feasible, specs)
        fallback = feasible[torch.argmin(penalties)]

        interior = raw_points[raw_means >= threshold + CONSTRAINT_START_MARGIN]
        interior_penalties = prior_information_penalty(interior, specs)
        restart_count = min(RECOMMENDATION_RESTARTS, len(interior))
        starts = interior[
            torch.topk(-interior_penalties, restart_count).indices
        ].unsqueeze(1)

    def accuracy_constraint(point: torch.Tensor) -> torch.Tensor:
        return model.posterior(point.unsqueeze(0)).mean.squeeze() - threshold

    try:
        candidate, _ = optimize_acqf(
            NegativePriorInformation(model, specs), bounds=unit_bounds(len(specs)), q=1,
            num_restarts=restart_count, raw_samples=None,
            batch_initial_conditions=starts,
            nonlinear_inequality_constraints=[(accuracy_constraint, True)],
            options={"batch_limit": 1, "maxiter": 200},
        )
    except CandidateGenerationError:
        print(
            "constrained recommendation optimization failed; "
            "using the best feasible Sobol candidate",
            flush=True,
        )
        return fallback

    candidate = candidate.squeeze(0)
    with torch.no_grad():
        candidate_mean = model.posterior(candidate.unsqueeze(0)).mean.squeeze()
        candidate_penalty = prior_information_penalty(candidate.unsqueeze(0), specs).squeeze()
        fallback_penalty = prior_information_penalty(fallback.unsqueeze(0), specs).squeeze()
        use_candidate = candidate_mean >= threshold and candidate_penalty <= fallback_penalty
    return candidate if bool(use_candidate) else fallback


def select_candidates(
    model: SingleTaskGP, specs: tuple[ParameterSpec, ...], mode: str, seed: int,
) -> tuple[dict, dict]:
    optimum = optimize_posterior_mean(model, len(specs))
    less_informative = optimize_less_informative(model, specs, optimum, seed)
    return (
        describe_candidate(model, optimum, specs, mode),
        describe_candidate(model, less_informative, specs, mode),
    )


def ard_summary(
    model: SingleTaskGP, specs: tuple[ParameterSpec, ...]
) -> tuple[dict[str, float], dict[str, float]]:
    values = model.covar_module.lengthscale.detach().reshape(-1).tolist()
    lengthscales = {spec.name: value for spec, value in zip(specs, values)}
    inverse = {name: 1.0 / value for name, value in lengthscales.items()}
    total = sum(inverse.values())
    normalized_inverse = {name: value / total for name, value in inverse.items()}
    return lengthscales, normalized_inverse


def summarize_records(
    records: list[dict], specs: tuple[ParameterSpec, ...], mode: str,
) -> dict:
    model = fit_surrogate(records)
    optimum, less_informative = select_evaluated_candidates(model, records, specs, mode)
    current_default = describe_candidate(
        model, torch.tensor(records[0]["x"], dtype=torch.double), specs, mode
    )
    current_default["evaluation_index"] = 0
    lengthscales, normalized_inverse = ard_summary(model, specs)
    noise_variance = observation_noise_variance(model)
    return {
        "current_default": current_default,
        "optimum": optimum,
        "less_informative": less_informative,
        "observation_noise_variance": noise_variance,
        "observation_noise_sd": math.sqrt(noise_variance),
        "ard_lengthscales": lengthscales,
        "normalized_inverse_lengthscales": normalized_inverse,
    }


def tune_mode(
    mode: str, seq_file: Path, dataset_index: int, args: argparse.Namespace,
) -> dict:
    """Run coverage, noisy BO, and adaptive recommendation for one mode."""
    specs = parameter_specs(mode)
    mode_index = args.modes.index(mode)
    seed = args.seed + 100_000 * dataset_index + 1_000 * mode_index

    current_default = evaluate_point(
        default_unit_point(specs), specs, mode, seq_file,
        args.mask, args.val, "current_default", 0,
    )
    records = [current_default]

    sobol = torch.quasirandom.SobolEngine(len(specs), scramble=True, seed=seed)
    for iteration, point in enumerate(sobol.draw(args.initial_points).double()):
        records.append(
            evaluate_point(
                point.tolist(), specs, mode, seq_file,
                args.mask, args.val, "sobol", iteration,
            )
        )

    for iteration in range(args.iterations):
        records.append(
            evaluate_point(
                propose_next_point(records), specs, mode, seq_file,
                args.mask, args.val, "bo", iteration,
            )
        )

    recommendation_history = []
    for iteration in range(args.recommendation_rounds):
        model = fit_surrogate(records)
        optimum, less_informative = select_candidates(
            model, specs, mode, seed + 10_000 + iteration
        )
        optimum_record = evaluate_point(
            optimum["x"], specs, mode, seq_file,
            args.mask, args.val, "recommend_optimum", iteration,
        )
        records.append(optimum_record)
        optimum["evaluation_index"] = len(records) - 1

        less_informative_record = evaluate_point(
            less_informative["x"], specs, mode, seq_file,
            args.mask, args.val, "recommend_less_informative", iteration,
        )
        records.append(less_informative_record)
        less_informative["evaluation_index"] = len(records) - 1
        recommendation_history.append(
            {"round": iteration, "optimum": optimum,
             "less_informative": less_informative}
        )

    return {
        "parameter_specs": [asdict(spec) for spec in specs],
        "records": records,
        "recommendation_history": recommendation_history,
        "summary": summarize_records(records, specs, mode),
    }


def resolve_seq_files(args: argparse.Namespace) -> list[Path]:
    if args.seq_files:
        files = args.seq_files
    elif args.seq_dir:
        files = sorted(args.seq_dir.glob("haps*.txt*"), key=get_seq_sort_key)
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
    parser.add_argument(
        "--initial_points",
        type=int,
        default=32,
        help="number of distinct initial Sobol evaluations (default: 32)",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=80,
        help="number of sequential qLogNEI evaluations (default: 80)",
    )
    parser.add_argument(
        "--recommendation_rounds",
        type=int,
        default=5,
        help="adaptive optimum/less-informative evaluation rounds (default: 5)",
    )
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
    if args.iterations < 0 or args.recommendation_rounds < 0:
        raise ValueError("iterations and recommendation_rounds must be nonnegative")


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
            "recommendation_rounds": args.recommendation_rounds,
            "seed": args.seed,
            "acquisition": {
                "name": "qLogNoisyExpectedImprovement",
                "mc_samples": ACQUISITION_MC_SAMPLES,
                "restarts": ACQUISITION_RESTARTS,
                "raw_samples": ACQUISITION_RAW_SAMPLES,
            },
        },
        "datasets": [],
    }
    for dataset_index, seq_file in enumerate(seq_files):
        dataset = {"seq_file": str(seq_file), "modes": {}}
        for mode in args.modes:
            dataset["modes"][mode] = tune_mode(mode, seq_file, dataset_index, args)
        experiment["datasets"].append(dataset)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(experiment, indent=0, allow_nan=False) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
