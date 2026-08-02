# DFCP Project Guide

This file is the working guide for future changes. Read it before editing the
project, then consult `notes.typ` and the relevant source files for details.
The C++ implementation is authoritative when the notes and code differ.

## Project Summary

DFCP stands for Discrete Fragmentation Coagulation Processes. The program fits
a sequence of related partitions to haplotype sequences across genomic loci.
At locus `l`, `R_l` partitions the sequences. Between adjacent loci, `R_l` is
fragmented into `Q_l`, then `Q_l` is coagulated into `R_(l+1)`.

The active implementation is the C++20 executable in `src/` and `include/`.
The files in `py/` are the deprecated prototype. They can help explain the
history of an algorithm, but new behavior belongs in C++ and should not be
ported from Python without checking it against the current model.

`notes.typ` is the main mathematical notebook. It includes derivations,
experiments, corrections to earlier derivations, and future ideas. In
particular, note the explicit corrections and warnings near the slice-sampling,
`d_l`, and maximization sections. It is not a line-by-line specification of the
current executable.

## Repository Map

### Root files

- `notes.typ`: model derivations, inference notes, metrics, and research TODOs.
- `CMakeLists.txt`: builds the single `dfcp` C++20 executable with OpenMP.
- `build.sh`: configures `build/`, links `compile_commands.json`, and builds.
- `README.md`: short project title only; use this guide for operational docs.
- `.clangd`: strict unused-include and missing-include diagnostics.
- `gp.ipynb`: separate Gaussian-process/Bayesian-optimization exploration. It
  is not part of the executable and contains large embedded notebook output.
- `data/`, `output/`, `build/`, `typstbuild/`: ignored inputs or generated
  artifacts, not source code.

### C++ headers

- `include/hyperparams.hpp`: dimensions and prior hyperparameters.
- `include/params.hpp`: moments for the continuous variational approximation.
- `include/clusters.hpp`: `R`/`Q` graph nodes, ownership, assignments, and
  cluster mutation interface.
- `include/max.hpp`: sequence reassignment and insertion entry points.
- `include/expect.hpp`: continuous-parameter update entry point.
- `include/elbo.hpp`: approximate ELBO entry point.
- `include/math.hpp`: second-order delta-method helpers.
- `include/tree.hpp`: reference-tree parsing and tree-based metrics.
- `include/util.hpp`: flat indexing, early stopping, held-out observations,
  modes, and numeric argument parsing.
- `include/json.hpp`: small write-only JSON builder used for stdout results.

### C++ sources

- `src/me.cpp`: executable orchestration, input and argument parsing, fitting,
  imputation, evaluation, timings, and JSON output.
- `src/clusters.cpp`: cluster creation/deletion, graph links, assignments,
  hard-emission indexes, and soft-emission counts.
- `src/max.cpp`: hard and soft sequencewise Viterbi maximization.
- `src/expect.cpp`: Laplace updates for `alpha`, `gamma_l`, and `d_l`.
- `src/elbo.cpp`: approximate ELBO and variational entropy.
- `src/math.cpp`: reusable delta approximations.
- `src/tree.cpp`: fastsimcoal tree parsing, Fitch parsimony, clade metrics,
  coordinate alignment, and DOT output.
- `src/util.cpp`: locus modes and numeric CLI parsers.

### Scripts

- `scripts/tune.py`: random hyperparameter search using soft-mode imputation.
- `scripts/parsimony.py`: evaluates simulated error files against true trees
  and plots excess parsimony.
- `scripts/viz.py`: reads one result JSON document from stdin and plots the
  training trace, parameter means, and timings.

### Deprecated Python

- `py/generate.py`: historical simulator and a concise executable description
  of the original generative process.
- `py/me.py`: original hard-emission Maximization-Expectation prototype.
- `py/impute.py`: old elementwise-masking imputation benchmark.

Do not extend `py/`. The C++ version has different storage, sequential
initialization, soft emissions, validation behavior, tree evaluation, and JSON
output.

## Mathematical Model

Let `N` be the number of haplotypes, `L` the number of loci, and `K` the number
of allele categories.

The partition process is:

1. `R_0 ~ CRP(N, alpha)`.
2. `Q_l ~ Frag(R_l, 0, d_l)` for `l = 0, ..., L-2`.
3. `R_(l+1) ~ Coag(Q_l, alpha / d_l, 0)`.

The priors used by the current executable are:

- `alpha ~ Gamma(tau_1, tau_2)` with shape and rate.
- `d_l ~ Beta(v_1, v_2)`.
- `gamma_l ~ Gamma(phi_1, phi_2)` with shape and rate.

Defaults are `tau=(1,1)`, `v=(1,1)`, and `phi=(2,2)`.

### Hard emissions

Each `R` cluster has one fixed allele. Existing clusters are feasible only
when their emission matches the observation. Missing observations do not
constrain the cluster. New-cluster emission probabilities integrate the
locus-level categorical probabilities under the symmetric Dirichlet prior.

### Soft emissions

Each `R` cluster has its own categorical distribution with symmetric
Dirichlet concentration `gamma_l`. The distribution is integrated out, so the
predictive likelihood for observed allele `k` in cluster `a` is based on

```text
(gamma_l + n_(a,l,k)) / (K gamma_l + n_(a,l,observed)).
```

A new soft cluster predicts each category with probability `1/K`. Missing
observations do not contribute to `nk` or `n_obs`.

### Inference

The executable alternates:

1. Maximization: remove one sequence at a time and put it back along the best
   complete `R-Q-R-...` path found by right-to-left Viterbi messages.
2. Expectation: update independent approximations for `alpha`, every `d_l`,
   and every `gamma_l`.
3. Evaluation: calculate the approximate ELBO and update early stopping.

Positive parameters use log-space Laplace approximations. Discounts use a
logit-space Laplace approximation. Several expectations use second-order delta
approximations. This is a MAP partition path with approximate continuous
posteriors, not a full posterior over partitions.

Early stopping has patience 2 and improvement tolerance `1e-3`. There is no
hard iteration limit.

## Core Data Structures And Invariants

The cluster graph is

```text
R_0 -> Q_0 -> R_1 -> Q_1 -> ... -> R_(L-1)
```

`Clusters::all_clusters` owns every node through `unique_ptr`. All pointers in
assignments, adjacency lists, and indexes are non-owning stable pointers.

Important invariants:

- Every active sequence belongs to exactly one `R` cluster at every locus.
- Every active sequence belongs to exactly one `Q` cluster at every transition.
- A `Q_l` is a fragment of one `R_l` and has exactly one child `R_(l+1)`.
- An `R` cluster may have multiple child fragments and multiple parent
  fragments.
- Empty clusters are detached and deleted immediately.
- `Cluster::n` counts all assigned sequences.
- In soft mode, `n_obs` excludes missing values and `nk[k]` counts observed
  allele `k` values only.
- In hard mode, only `R` clusters have emissions; `Q` emissions are `-1`.
- `r_assign` is flat `[N][L]`; `q_assign` is flat `[N][L-1]`.
- Use `idx2d` for flat indexing rather than reproducing index arithmetic.

`Clusters` stores a reference to `HyperParams`. `HP.N` is deliberately changed
during sequential initialization and validation insertion, so understand the
current value before using it.

## Runtime Data Flow

`src/me.cpp` performs these steps:

1. Read one haplotype per line into a flat row-major `int8_t` vector.
2. Infer `N`, `L`, and `K` from the input.
3. Parse all optional arguments as option/value pairs.
4. Optionally split whole sequences into training and validation groups, then
   mask alleles only within held-out sequences.
5. Initialize parameters and clusters, either as one block or by adding
   sequences through Viterbi.
6. Run Maximization-Expectation until early stopping.
7. Optionally insert held-out sequences and score imputation.
8. Restore masked values for evaluation.
9. Optionally align loci to reference trees and compute tree metrics.
10. Compute adjacent-locus co-clustering IoU, mean cluster count, and purity.
11. Write diagnostics to stderr and one JSON object to stdout.

Validation reorders sequences internally. `raw_to_split_idxs` maps each leaf
index in the original input/tree order to its current row in `x` and
`r_assign`. Every reference-tree metric must apply this map.

## Input Formats

### Haplotype sequence file

- One sequence per line.
- One ASCII digit per allele.
- All lines are expected to have the same nonzero length.
- Alleles are limited to categories `0` through `9`.
- `K` is `max_observed_digit + 1`; labels need not be contiguous.
- The executable reads plain text. A filename containing `.gz` does not imply
  decompression and the supplied examples with such names are plain text.
- Internal missing values use `-1`, but `-1` is not accepted in input files.

### Variant position file

The parser expects an integer count followed by comma-separated integer
positions. `--variant_start_pos` is a zero-based index into this list, not a
genomic coordinate, despite the option name.

### Tree file

The parser targets fastsimcoal/NEXUS marginal-tree output. It skips three
header lines, reads records containing `pos_`, expects binary trees with the
spacing used by the example data, converts one-based leaf labels to zero-based
indices, and ignores branch lengths.

Trees are represented as

```cpp
std::unordered_map<int, std::pair<int, int>>
```

Internal node IDs are negative, the root is `-1`, and leaves are nonnegative
sequence indexes. A node absent from the map is treated as a leaf.

At a recombination boundary, `get_tree_idxs` selects the new tree for a variant
at exactly that position.

## Command-Line Interface

Build and run from the repository root:

```bash
./build.sh
./build/dfcp SEQUENCE_FILE [OPTION VALUE]...
```

Every option requires a value, including booleans. There is no `--help` path.

- `--tau_1`, `--tau_2`: Gamma shape/rate for `alpha`.
- `--v_1`, `--v_2`: Beta shapes for each `d_l`.
- `--phi_1`, `--phi_2`: Gamma shape/rate for each `gamma_l`.
- `--val`: probability of holding out a whole sequence.
- `--mask`: probability of masking each allele in held-out sequences.
- `--tree`: fastsimcoal marginal-tree file.
- `--variant_pos`: variant-position file.
- `--variant_start_pos`: starting index in the position list.
- `--tree_vis`: DOT output path for the first eight loci.
- `--clade_beta`: symmetric Beta shape for cluster-to-clade importance
  weights; defaults to `1` and must be at least `1`.
- `--soft`: enabled only when the value is exactly `1`.
- `--block_init`: enabled only when the value is exactly `1`.

A representative hard-mode run is:

```bash
./build.sh && ./build/dfcp \
  data/examples/simulated/SIM1_LEN500_NHAPS100/haps_SIMOUT_1.txt.gz_SIMOUT_14572-15071.txt \
  --tree data/examples/simulated/SIM1_LEN500_NHAPS100/ex_0_pop_1_1_true_trees.trees_1_14572_500.trees \
  --variant_pos data/examples/simulated/SIM1_LEN500_NHAPS100/variant_pos.txt \
  --variant_start_pos 14572 \
  --mask 0.2 --val 0.2
```

Validation is randomized through `std::random_device`; repeated runs need not
produce identical partitions or metrics.

## Output

The executable emits progress and human-readable metrics to stderr. Stdout is
one JSON object suitable for scripts.

Always-present fields include:

- `seq_file`
- `train_log`
- `params.mu_alpha`, `params.mu_gamma`, `params.mu_d`
- `mean_iou`, `mean_emission_iou`, `t_iou`
- `mean_clusters`
- `cluster_purity`

Validation adds held-out counts, DFCP and mode imputation accuracy, and timing.
Tree evaluation adds mean excess parsimony, emission excess parsimony,
`clade_beta`, importance-weighted `mean_clade_iou`, `t_parsimony`, and
`t_clade_iou`.

Cluster assignments are not currently serialized.

## Evaluation Metrics

### Excess parsimony

Fitch's algorithm assigns possible labels to internal nodes while minimizing
edge changes. A partition into `k` clades has minimum score `k-1`, so the
reported excess is

```text
Fitch score - (k - 1).
```

Zero means every inferred cluster can be represented as a tree clade.

### Importance-weighted cluster-to-clade IoU

For inferred cluster `C` and every clade `S` represented by a tree node,

```text
J(S, C) = |S intersect C| / |S union C|.
```

The per-cluster score is `max_S J(S,C)`. Leaves and the root are included as
singleton and full-tree clades. Although `notes.typ` calls this a Jaccard
distance, the displayed formula and the implementation are similarities:
higher is better. A distance would be one minus the reported value.

For cluster size `c`, define `z=(c-1)/(N-1)`. Given the symmetric shape
`clade_beta`, the aggregate metric uses weights proportional to:

```text
[z (1-z)]^(clade_beta - 1)
```

Beta normalizing constants and the interval scaling factor cancel between the
weighted numerator and denominator. The aggregate is over all clusters at all
loci, as written in `notes.typ`, rather than a mean of per-locus means.
Shape `1` gives every cluster equal weight. Shapes greater than `1` increasingly
downweight singleton and full-sample clusters. Shapes below `1` are rejected
because their density is infinite at those endpoints. If every cluster has
zero weight, the executable reports `1`, since those endpoint clusters match
leaf or root clades exactly.

The clade maximization is computed for all clusters in one postorder traversal
per locus. Each subtree carries counts by cluster label and child maps are
merged small-to-large. A leaf initializes the singleton-clade score. At an
internal node, a cluster is evaluated only if it occurs in both children. If
it occurs in one child only, the parent has the same intersection and extra
non-cluster leaves, so its IoU cannot beat that child. This avoids explicitly
evaluating every cluster at every clade while remaining exact.

### Adjacent-locus IoU

`mean_iou` compares the relations "sequence pair is co-clustered" at adjacent
loci. It is not a mean of Jaccard scores between individual clusters.
`mean_emission_iou` applies the same pairwise relation to equal observed
emissions.

### Purity

Hard clusters are emission-pure by construction, so hard-mode purity is `1`.
Soft-mode purity sums each cluster's majority observed allele count and divides
by `L * n_train_seqs`.

## Build And Dependencies

Required for the executable:

- CMake 3.20 or newer.
- A C++20 compiler with `<format>` support.
- OpenMP.
- Boost headers for special functions, logistic sigmoid, and Brent minimization.

The build enables `-Wall -Wextra -Wpedantic -O3`. There is currently no
install target, library target, automated test target, or CI configuration.

The active plotting scripts require Python and Matplotlib. Historical and
notebook work may additionally require NumPy, SciPy, Graphviz, Jupyter, and
line-profiler. There is no tracked Python package or lock manifest.

## Testing And Verification

After a C++ change, at minimum run:

```bash
./build.sh
```

For changes involving fitting or tree metrics, run the representative command
above. Redirect stdout only when JSON contents are not under test; diagnostics
remain visible on stderr.

For tree metric changes, use small hand-built trees or a temporary standalone
test outside the repository and compare against brute-force enumeration. Useful
cases are:

- An exact two-leaf clade, expected maximum IoU `1`.
- A cluster split across two sibling clades.
- A singleton cluster, expected maximum IoU `1` at its leaf.
- The full sample cluster, expected maximum IoU `1` at the root.
- An internal ancestor with members in only one child, which must be dominated
  by that child.
- Random small binary trees and partitions checked against every explicit
  descendant-leaf set.
- Validation-enabled runs, to exercise `raw_to_split_idxs`.

There is no deterministic seed option, so validation runs are smoke tests
unless the input or executable is extended to control the RNG.

## Known Caveats

- Empty input and ragged input lines are not explicitly validated before all
  indexing operations.
- The executable assumes `L >= 2`; adjacent-locus metrics divide by `L-1`.
- Numeric argument parsing accepts a valid numeric prefix followed by junk.
- Validation fractions are not range-checked.
- The transformed Laplace searches use fixed bounds `[-10,10]` and fixed
  precision.
- Unordered containers and exact score ties can make fitting nondeterministic.
- Tree parsing is intentionally specialized and does not validate arbitrary
  Newick trees.
- DOT output is written as an eight-locus graph. Inputs shorter than eight loci
  do not reach the expected closing-locus condition.
- Soft purity after validation insertion includes inserted clusters/counts in
  the numerator but uses `n_train_seqs` in the denominator.
- `Json` rejects NaN and infinity, so metrics must define finite edge-case
  behavior.
- `notes.typ` contains both current soft-model work and superseded hard-model
  derivations. Check nearby correction notes and the C++ implementation.

## C++ Style Guide

Match the existing code unless a local cleanup is required for correctness.

- Use C++20 and standard-library facilities directly.
- Keep public declarations in the matching header and implementation in the
  matching `.cpp` file.
- Prefer small free functions and simple structs. Do not add class hierarchies,
  templates, or abstractions unless they remove real repeated complexity.
- Keep an operation in one function when it is used once and remains easy to
  read. Extract a helper when recursion, reuse, or focused testing makes the
  main path substantially clearer.
- Use `unique_ptr` for ownership and raw pointers only for explicit non-owning
  graph links, matching `Clusters`.
- Prefer flat vectors and `idx2d` in hot sequence/locus code.
- Use references for required non-owning inputs and `const` wherever practical.
- Follow the mathematical names already used by the project: `a` for `R`
  clusters, `b` for `Q` clusters, `l` for locus, and `i` for sequence.
- Use descriptive names for metric accumulators and data that are not part of
  the mathematical derivation.
- Keep control flow direct. Prefer early validation and early returns over deep
  nesting.
- Braces are required for multi-line blocks. Existing one-line throw guards are
  acceptable.
- Use four spaces for indentation and keep blank-line spacing consistent with
  neighboring code.
- Include the standard header that declares every directly used symbol. Keep
  standard headers before project headers and satisfy strict clangd include
  diagnostics.
- Do not introduce a project namespace solely for style consistency; the
  current code uses global declarations.
- Comments should explain an invariant, derivation, pruning argument, or
  non-obvious representation. Do not narrate straightforward assignments.
- Throw `std::invalid_argument` for invalid caller/input contracts and
  `std::runtime_error` for broken runtime or model invariants, following the
  existing distinction.
- Preserve stderr for diagnostics and stdout for the single JSON result.
- Add JSON fields with clear semantic names. Keep existing fields stable unless
  a deliberate interface change is required.
- Avoid speculative compatibility layers. Update the active interface directly
  when no external compatibility requirement is known.
- Optimize only where the data shape justifies it, but preserve obvious
  asymptotic improvements in tree and pairwise code. Document exact pruning
  arguments close to the implementation.

## Python Script Style

Active scripts are small command-line programs rather than a package.

- Use `argparse` for user inputs and `subprocess.run(..., check=True)` when a
  failed command should stop the script.
- Use `pathlib.Path` for filesystem traversal.
- Parse executable stdout as JSON; do not scrape stderr diagnostics.
- Keep scripts runnable from the repository root, matching current paths.
- Do not add dependencies for tasks the standard library handles clearly.
- Do not add new logic to `py/`; place maintained utilities in `scripts/`.

## Change Checklist

Before editing:

- Identify whether the behavior is hard mode, soft mode, or shared.
- Trace sequence indexing through validation and `raw_to_split_idxs`.
- Check graph ownership and deletion if cluster pointers are involved.
- Check `notes.typ` for a correction near the relevant derivation.

Before finishing:

- Build with `./build.sh` and resolve all warnings.
- Run an end-to-end command that exercises the changed path.
- Check stdout remains valid JSON and stderr contains diagnostics only.
- Test edge cases of any new numerical denominator or empty collection.
- Update this guide when architecture, interfaces, formulas, or caveats change.
