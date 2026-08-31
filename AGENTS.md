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
- `CMakeLists.txt`: builds the `impute` and `eval_impute` C++20 executables.
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
- `include/seq_array.hpp`: sequence-major bitpacked observations and the binary
  sequence-file reader.
- `include/io.hpp`: shared little-endian binary I/O and atomic output writing.
- `include/impute_io.hpp`: streamed fixed-point imputation probability I/O.
- `include/clusters.hpp`: `R`/`Q` graph nodes, ownership, assignments, and
  cluster mutation interface, including the emission-mode enum.
- `include/max.hpp`: sequence reassignment and insertion entry points.
- `include/expect.hpp`: continuous-parameter update entry point.
- `include/elbo.hpp`: approximate ELBO entry point.
- `include/math.hpp`: second-order delta-method helpers.
- `include/tree.hpp`: reference-tree parsing and tree-based metrics.
- `include/util.hpp`: flat indexing, early stopping, held-out observations,
  modes, and numeric argument parsing.
- `include/json.hpp`: small write-only JSON builder used for stdout results.

### C++ sources

- `src/impute.cpp`: executable orchestration, input and argument parsing, fitting,
  imputation, timings, and JSON output.
- `src/seq_array.cpp`: binary sequence loading and the 64-by-64 bit transpose
  from locus-major file words to sequence-major memory words.
- `src/io.cpp`: exact binary reads/writes, endian conversion, and atomic file output.
- `src/impute_io.cpp`: little-endian imputation probability I/O.
- `src/eval_impute.cpp`: pooled minor-allele r-squared and accuracy across
  materialized VCF windows.
- `src/clusters.cpp`: cluster creation/deletion, graph links, assignments,
  hard-emission indexes, and soft-emission counts.
- `src/max.cpp`: hard and soft sequencewise Viterbi maximization.
- `src/expect.cpp`: Laplace updates for `alpha`, `gamma_l`, and `d_l`.
- `src/elbo.cpp`: approximate ELBO and variational entropy.
- `src/math.cpp`: reusable delta approximations.
- `src/tree.cpp`: fastsimcoal tree parsing, Fitch parsimony, clade metrics,
  coordinate alignment, and DOT output.
- `src/util.cpp`: emission counts and locus modes.

### Scripts

- `scripts/init.py`: compares initialization methods across DFCP modes, training
  phases, masks, and simulated datasets, saves reusable JSON results, and can
  rebuild its interactive Plotly report without rerunning DFCP.
- `scripts/dfcp.py`: shared Python interface for building and running the C++
  executable and parsing its JSON output.
- `scripts/tune.py`: JSON-only noisy BoTorch optimization of interpretable prior
  parameters for hard, noisy, and soft PBWT-initialized DFCP, independently for
  each sequence file, with learned observation noise and adaptive recommendations.
- `scripts/tune_viz.py`: reads tuning JSON and builds the interactive Plotly
  hyperparameter-sensitivity report without rerunning DFCP.
- `scripts/seq_file_name.py`: parses injected bit-flip and switch-error rates from
  simulated haplotype filenames for consistent labels and ordering.
- `scripts/plotly_html.py`: writes one shared `docs/assets/plotly.min.js` for generated
  reports and can externalize the Plotly bundle from an existing report.
- `scripts/parsimony.py`: evaluates simulated error files against true trees
  and plots excess parsimony.
- `scripts/viz.py`: reads one result JSON document from stdin or a file and
  builds an interactive Plotly report; it renders DFCP's DOT tree output as
  locus-selectable Graphviz SVGs with zoom controls.
- `scripts/1000g_phase3_v5b/window.py`: writes aligned, equal-locus overlapping
  reference and target VCF windows before DFCP bitpacking.
- `scripts/1000g_phase3_v5b/prep_data.sh`: runs the complete 1000 Genomes data
  preparation pipeline.
- `scripts/1000g_phase3_v5b/impute.sh`: builds and imputes every materialized
  window, writing `probs.bin` and `impute.json` in each directory.
- `scripts/1000g_phase3_v5b/window_viz.py`: interactive physical/genetic
  window-selection report using the same boundary formula as `window.py`.
- `scripts/impute_viz.py`: plots pooled imputation r-squared and accuracy by
  reference minor-allele count from `eval_impute` TSV output.

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

1. Expectation: update independent approximations for `alpha`, every `d_l`,
   and every `gamma_l`.
2. Maximization: remove one sequence at a time and put it back along the best
   complete `R-Q-R-...` path found by right-to-left Viterbi messages.
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

`Clusters::all_clusters` is an ID-indexed vector that owns every node through
`unique_ptr`; deleted IDs leave reusable null slots. All pointers in assignments,
adjacency lists, and indexes are non-owning stable pointers. Every cluster has a
`uint32_t` ID that is stable for its lifetime and returned to a free list after
deletion. Viterbi and forward-backward R messages use dense vectors indexed by
this ID; proposed-new-cluster messages remain separate by locus. `ViterbiBuffers`
owns the reusable Viterbi messages and path. `FwdBkwdBuffers` owns reusable
forward and backward message arrays. Imputation also overwrites one probability
row for every target sequence in both inference modes.
Active `R` and `Q` clusters and the hard-emission index are dense pointer
vectors. Empty-cluster deletion finds the pointer linearly and uses
swap-and-pop; the vectors average only a few entries per locus.

Important invariants:

- Every active sequence belongs to exactly one `R` cluster at every locus.
- Every active sequence belongs to exactly one `Q` cluster at every transition.
- A `Q_l` is a fragment of one `R_l` and has exactly one child `R_(l+1)`.
- Q clusters store their parent and child directly; only R clusters use the
  `parents` and `children` vectors. A Q link may briefly be null when its adjacent
  empty R cluster is deleted first during sequence removal.
- An `R` cluster may have multiple child fragments and multiple parent
  fragments.
- Empty clusters are detached and deleted immediately.
- A live cluster ID is unique. Reused IDs are assigned only after the previous
  cluster is detached and deleted.
- Inference reads messages only for current candidates, whose ID-indexed entries
  were overwritten during the current pass.
- Use `Clusters::get_matching_as` for R candidates: hard observed emissions
  select `rs_by_emit`, while missing, noisy, and soft emissions select all `rs`.
- `Cluster::n` counts all assigned sequences.
- In soft mode, `n_obs` excludes missing values and `nk[k]` counts observed
  allele `k` values only.
- In hard mode, only `R` clusters have emissions; `Q` emissions are `-1`.
- `SeqArray::x` is bitpacked sequence-major `[N][ceil(L/64)]`.
- `r_assign` is flat `[N][L]`; `q_assign` is flat `[N][L-1]`.
- Use `idx2d` for flat indexing rather than reproducing index arithmetic.

`Clusters` stores a reference to `HyperParams`. `HP.N` is deliberately changed
during sequential initialization and validation splitting, so understand the
current value before using it.

## Runtime Data Flow

`src/impute.cpp` performs these steps:

1. Read locus-major bitpacked reference and target files and transpose each
   into an in-memory sequence-major `SeqArray`.
2. Read the map from full reference loci to compact observed-target columns.
3. Parse all optional arguments as option/value pairs.
4. Initialize parameters and clusters, either as one block, with PBWT groups,
   or by adding sequences through Viterbi.
5. Unless `--init_only 1` is set, run Expectation-Maximization until early
   stopping.
6. Impute target loci absent from the observed-target map without inserting
   target sequences into the fitted reference cluster graph, streaming one
   fixed-point probability row at a time.
7. Write diagnostics to stderr and one JSON object to stdout.

## Input Formats

### Haplotype sequence file

- Sequence files are binary and begin with the four ASCII magic bytes `DFCP`,
  equivalent to the little-endian 32-bit integer `0x50434644`.
- The magic is followed by little-endian 32-bit `N` and `L` values.
- The payload is locus-major `[L][ceil(N/64)]` little-endian 64-bit words.
  Allele `(i,l)` is bit `i%64` of payload word
  `l*ceil(N/64)+i/64`, with the least-significant bit first.
- Every locus begins at a 64-bit boundary. Unused high bits in its final word
  must be zero, and trailing bytes are rejected.
- The reader transposes 64-by-64 bit tiles into sequence-major
  `[N][ceil(L/64)]` `SeqArray` storage. Every in-memory sequence therefore also
  begins at a 64-bit boundary.
- The file format supports only alleles `0` and `1`; missing values are not
  represented. `N` and `L` must both be nonzero and fit in a positive `int`.

### Observed loci file

- The file contains one zero-based reference-locus index per observed target
  column, in the same order as columns in the target sequence file.
- It must contain exactly `target.L` unique indexes, each in `[0, reference.L)`.
- `dfcp_prep.sh` creates it for each VCF window by matching `CHROM`, `POS`,
  `REF`, and `ALT` and rejects missing, duplicate, or reordered target variants.

### Imputation probability file

- The four-byte magic is `DFIP`, followed by little-endian `uint32` target
  sequence count `N` and masked-locus count `M`.
- The payload is row-major `[N][M]` little-endian `uint16` values encoding
  `round(P(allele 1) * 65535)`.
- The masked-locus order is the complement of the observed-loci file in
  increasing reference-locus order.
- The writer streams one target row and atomically replaces the output after
  all rows succeed.
- Allele 1 is not guaranteed to be the reference-panel minor allele. Compute
  minor count as `min(allele_1_count, N - allele_1_count)`.

### Imputation evaluation TSV

- `eval_impute` reads a window root containing `windows.tsv` and evaluates each
  generated window that contains `probs.bin`.
- Each global locus is retained from the available window in which it is
  farthest from an edge. Ties go to the lower window index.
- Reference `AC` and `AN` are read from each window's `ref.vcf.gz` with
  `bcftools`. Probability and truth alleles are flipped when REF is minor.
- The output columns are `mac`, `n_loci`, `n_predictions`, `r2`, and `accuracy`.
  Statistics are pooled over all retained target alleles in each exact MAC bin.
- An r-squared value of `-1` means the truth or probability had zero variance
  within that MAC bin.

### Variant position file

This format is used by the currently inactive evaluation sources, not by
`dfcp`'s active command line.

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
std::unordered_map<int, CoalNode>
```

Internal node IDs are negative, the root is `-1`, and leaves are nonnegative
sequence indexes. A node absent from the map is treated as a leaf. `CoalNode`
stores both child indexes and their branch lengths.

At a recombination boundary, `get_tree_idxs` selects the new tree for a variant
at exactly that position.

## Command-Line Interface

Build and run from the repository root:

```bash
./build.sh
./build/impute REF_FILE TARGET_FILE OBSERVED_LOCI_FILE PROB_FILE [OPTION VALUE]...
```

Every option requires a value, including booleans. There is no `--help` path.

- `--tau_1`, `--tau_2`: Gamma shape/rate for `alpha`.
- `--v_1`, `--v_2`: Beta shapes for each `d_l`.
- `--phi_1`, `--phi_2`: Gamma shape/rate for each `gamma_l`.
- `--mode`: emission model, one of `hard`, `noisy`, or `soft`; defaults to
  `hard`.
- `--lambda_1`, `--lambda_2`: Beta shapes for the noisy-emission error rate.
- `--init`: initialization method, one of `viterbi`, `block`, or `pbwt`;
  defaults to `pbwt`.
- `--pbwt_match_len`: PBWT initialization match length; defaults to `5`.
- `--pbwt_match_curr`: include the current locus in PBWT matching only when the
  value is exactly `1`; defaults to enabled.
- `--init_only`: skip ME training only when the value is exactly `1`.
- `--viterbi_impute`: use the Viterbi path rather than forward-backward
  imputation only when the value is exactly `1`.

A prepared 1000 Genomes invocation is:

```bash
./build.sh && ./build/impute \
  data/1000g_phase3_v5b/windows/window_0000/ref.bin \
  data/1000g_phase3_v5b/windows/window_0000/target_observed.bin \
  data/1000g_phase3_v5b/windows/window_0000/observed_loci.txt \
  data/1000g_phase3_v5b/windows/window_0000/probs.bin \
  --init pbwt --viterbi_impute 1
```

Evaluate the probabilities with:

```bash
./build/eval_impute \
  data/1000g_phase3_v5b/windows \
  output/imputation.tsv
.venv/bin/python scripts/impute_viz.py \
  output/imputation.tsv \
  --output output/imputation.html
```

The evaluator expects each completed window's probability file to be named
`probs.bin`.

## Output

The executable emits progress and human-readable metrics to stderr. Stdout is
one JSON object suitable for scripts.

Always-present fields include:

- `ref_file`, `target_file`, and `observed_loci_file`
- `t_init` and `t_impute`

Unless `--init_only 1` is set, output also includes `train_log` and the fitted
parameter moments under `params`. Each training record includes `mean_nR`, the
mean number of active R clusters per locus. Probabilities are written to
`PROB_FILE` and are not accumulated in JSON. Cluster assignments are not serialized.

The cluster and tree metrics below describe intended behavior in inactive
evaluation sources. Imputation r-squared and accuracy are active in
`eval_impute`; they are not emitted by `dfcp`.

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
zero weight, the executable reports `-1` to mark the aggregate as undefined.

The clade maximization uses a postorder traversal for each cluster. With
validation enabled, held-out leaves contribute neither intersections nor clade
sizes, which is equivalent to evaluating descendant-leaf sets on the tree
induced by training leaves. Unary ancestors created by pruning are dominated
by their retained child.

### Adjacent-locus IoU

`mean_iou` compares the relations "sequence pair is co-clustered" at adjacent
loci. It is not a mean of Jaccard scores between individual clusters.
`mean_emission_iou` applies the same pairwise relation to equal observed
emissions.

### Purity

Hard clusters are emission-pure by construction, so hard-mode purity is `1`.
Soft-mode purity sums each cluster's majority observed allele count and divides
by `L * n_train_seqs`.

### Minor-allele imputation r²

For binary data, the minor allele at each masked locus is defined from the
reference panel. `eval_impute` pools predicted minor-allele probabilities and
true minor-allele indicators across every retained target allele at the same
reference MAC, then reports their squared Pearson correlation. It likewise
reports pooled hard-call accuracy rather than a mean of per-locus accuracies.
The online covariance calculation returns `-1` when either value is constant
within a MAC bin.

## Build And Dependencies

Required for the executable:

- CMake 3.20 or newer.
- A C++20 compiler with `<format>` support.
- OpenMP.
- Boost headers for special functions, logistic sigmoid, and Brent minimization.
- `bcftools` for windowed imputation evaluation.

The build enables `-Wall -Wextra -Wpedantic -O3`. There is no install target,
library target, automated test target, or CI configuration.

The active analysis scripts require Python, Matplotlib, Plotly, PyTorch,
BoTorch, and GPyTorch. Historical and notebook work may additionally require
NumPy, SciPy, Graphviz, Jupyter, and line-profiler. There is no tracked Python
package or lock manifest.

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
- Validation-enabled runs, to exercise `train_idxs` and pruning of held-out
  tree leaves.

Use a fixed nonzero `--seed` for deterministic validation and mask selection.

## Known Caveats

- The bitpacked sequence format does not represent missing or multiallelic
  observations. Missing target loci are represented by the separate mapping
  from full reference loci to columns in the compact target `SeqArray`.
- The executable assumes `L >= 2`; adjacent-locus metrics divide by `L-1`.
- Numeric argument parsing accepts a valid numeric prefix followed by junk.
- The transformed Laplace searches use fixed bounds `[-10,10]` and fixed
  precision.
- Unordered containers and exact score ties can make fitting nondeterministic.
- Tree parsing is intentionally specialized and does not validate arbitrary
  Newick trees.
- DOT output is written as a 16-locus graph. Inputs shorter than 16 loci
  do not reach the expected closing-locus condition.
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
- Establish pipeline invariants once where data is generated or first read. Do
  not repeat those checks downstream; only validate independently produced
  inputs when their compatibility is required for safe interpretation.
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
- Trace sequence indexing through validation and `train_idxs`.
- Check graph ownership and deletion if cluster pointers are involved.
- Check `notes.typ` for a correction near the relevant derivation.

Before finishing:

- Build with `./build.sh` and resolve all warnings.
- Run an end-to-end command that exercises the changed path.
- Check stdout remains valid JSON and stderr contains diagnostics only.
- Test edge cases of any new numerical denominator or empty collection.
- Update this guide when architecture, interfaces, formulas, or caveats change.
