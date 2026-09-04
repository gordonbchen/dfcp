# DFCP TODO

## Current objective

Separate model fitting from imputation, serialize frozen models and final R assignments, and restore
`eval_clusters` against a reproducible fastsimcoal fixture. Complete and test one stage before starting the
next.

## Current state

- `scripts/fsc_sim/run.sh` reproducibly runs fastsimcoal with an explicit seed, then
  `prep_data.py` converts its haploid `.gen` file into `ref.bin` and `variant_pos.txt`. The
  `.gen` file is the single source for both alleles and positions.
- The active `impute` executable still trains and imputes in one invocation. Models and cluster assignments
  are not serialized.
- `src/eval_clusters.cpp` is legacy code, uses removed interfaces, and is not built by CMake.
- The 1000 Genomes VCF preparation, windowed imputation, `DFIP` output, pooled `eval_impute`, and
  `impute_viz.py` pipeline are active.
- Viterbi and forward-backward messages, paths, and probability rows reuse dense buffers indexed by stable
  cluster IDs. Active cluster collections are dense vectors, while ID-indexed ownership uses reusable slots.

## Decisions

### fastsimcoal preparation

- Use the `.gen` file rather than parsing positions from `.arp` and alleles from `.gen`. The two position
  lists are identical in the current simulation, and one source cannot drift out of alignment with itself.
- Preserve repeated positions. The current fastsimcoal output has nine, and the tree lookup supports
  nondecreasing positions.
- Keep all haplotypes in their generated order so assignment row `i` corresponds to fastsimcoal tree leaf
  `i + 1` before the tree parser converts leaves to zero-based indexes.
- Do not create a synthetic VCF for cluster evaluation. Add one later only if a simulated end-to-end VCF
  imputation fixture becomes useful.

### Executable split

- `train` will take a reference sequence file and model output path. It will optionally write final R
  assignments when an assignment path is provided.
- `impute` will take a frozen model, compact target sequence file, observed-loci file, and probability output
  path. It will not fit or mutate the model.
- Keep training diagnostics and one compact JSON result on stdout. Imputation should report only its own
  inputs, dimensions, method, and timing.
- Update the 1000 Genomes scripts and maintained Python wrappers only after both executables work directly.

### Binary outputs

- `DFCM` is a frozen inference model. It stores dimensions, emission mode, fitted parameters, file-local
  dense cluster records, and R-Q graph links. It does not store training assignments.
- Derived containers such as per-locus active vectors and hard-emission indexes are rebuilt when loading.
- `DFCA` stores final R assignments as magic, `uint32 N`, `uint32 L`, then row-major `[N][L]` little-endian
  `uint32` cluster IDs.
- Do not store `next_cluster_id` in `DFCA`. It is an allocation high-water mark shared by R and Q nodes, not
  the number of live R clusters. `eval_clusters` will remap IDs densely at each locus.
- Use a path-valued training option such as `--output_clusters FILE`; omitting it disables assignment output.

## Current plan: restore cluster evaluation

### 1. Prepare fastsimcoal data

Status: complete.

- Run fastsimcoal with `-G`, `-T`, and an explicit seed. Suppress the unused Arlequin output;
  retain the haploid genotype table and true marginal trees.
- Convert `.gen` positions to the count-plus-comma-separated format consumed by the tree evaluator.
- Bitpack haploid allele columns locus-major into the existing `DFCP` sequence format.
- Reject incompatible `.gen` inputs at this I/O boundary: unexpected columns, multiple chromosomes,
  decreasing positions, nonbinary alleles, or empty data.
- Verify every decoded bit and position against `data/fsc/ex_0_pop_1/ex_0_pop_1_1_1.gen`.

### 2. Specify and implement frozen model I/O

Status: next.

- Write the exact `DFCM` byte layout beside its reader and writer before changing either executable.
- Store every fitted `Params` value used by Viterbi and forward-backward inference.
- Serialize only state needed for frozen inference: R and Q locus, size, emission/count state, and
  graph links, using dense file-local IDs. Each Q record can identify its one parent and child R;
  do not also serialize redundant R adjacency arrays.
- Load ownership first, resolve Q links and rebuild R adjacency second, then rebuild active locus
  vectors and `rs_by_emit`.
- Validate dimensions, IDs, cluster types, loci, and links once while reading.
- Prove that inference from the in-memory model and the saved-then-loaded model produces byte-identical
  probability output in hard, noisy, and soft modes.

### 3. Split training and imputation

Status: pending after model I/O.

- Build `train REF_FILE MODEL_FILE [OPTION VALUE]...` from the current initialization and ME loop.
- Build `impute MODEL_FILE TARGET_FILE OBSERVED_LOCI_FILE PROB_FILE [OPTION VALUE]...` around frozen model
  loading and the existing Viterbi or forward-backward target pass.
- Keep initialization, prior, batching, and stopping options on `train`; keep only inference options on
  `impute`.
- Remove fitted parameter arrays from training JSON once `DFCM` is authoritative.
- Update `scripts/1000g_phase3_v5b/impute.sh`, `scripts/dfcp.py`, `scripts/init.py`, and
  `scripts/tune.py` to the split interface where they are still maintained.

### 4. Add optional R-assignment output

Status: pending after `train` exists.

- Add shared little-endian `uint32` array I/O if the model and assignment formats both need it.
- Stream or write one sequence row at a time through `AtomicBinaryWriter`.
- Write assignments after the final fitted state, including after initialization when `--init_only 1` is
  used.
- Add a `DFCA` reader that checks its dimensions and payload size at the I/O boundary.
- Verify every written ID equals `clusters.r_assign[idx2d(i, l, L)]->id` before relying on the file in
  evaluation.

### 5. Rewrite `eval_clusters`

Status: pending after assignment output.

- Make it a standalone executable taking `ref.bin`, assignments, `variant_pos.txt`, and the fastsimcoal true
  tree file. The prepared position file already has exactly `L` entries, so remove `variant_start_pos`.
- Dense-remap raw assignment IDs independently at every locus and derive cluster counts and sizes from the
  assignments.
- Restore adjacent-locus cluster IoU, mean clusters per locus, excess parsimony, importance-weighted
  cluster-to-clade IoU, purity, emission IoU, emission parsimony, and emission clade IoU.
- Keep tree visualization optional and separate from scalar evaluation output.
- Replace avoidable `O(LN^2)` pair enumeration with contingency counts if that can be done without obscuring
  the metric.
- Emit one JSON document to stdout and diagnostics to stderr.

### 6. Verify the simulated pipeline

Status: pending.

- Train on all 100 generated haplotypes so assignment rows map directly to the 100 tree leaves.
- Evaluate all 13,624 generated variants against `ex_0_pop_1_1_true_trees.trees`.
- Test tree selection immediately before, at, and after recombination boundaries, including repeated variant
  positions.
- Compare tree routines with brute force on exact clades, sibling-split clusters, singleton and full-sample
  clusters, dominated ancestors, and random small binary trees.
- Record exact commands, dimensions, runtime, peak memory, and the resulting metrics.

## Later backlog

### Automated regression coverage

- Add a small VCF fixture covering overlapping windows, `--n-generate`, manifest rows, aligned records,
  observed-locus indexes, bitpacked alleles, and output dimensions.
- Add focused `DFIP` and pooled-evaluator tests for fixed-point decoding, overlap ownership, minor-allele
  orientation, mismatched dimensions, constant bins, and the `0.5` hard-call threshold.
- Keep one deterministic end-to-end fixture covering hard, noisy, and soft modes and both imputation methods.
- Verify masked indexes are the complement of observed indexes and VCF `min(AC, AN-AC)` agrees with decoded
  reference counts.

### Resumable training checkpoints

- Extend model serialization with an explicit assignments-present flag only after frozen loading works.
- Store both R and Q assignments and all mutable state needed to resume training exactly.
- Verify that an uninterrupted run and a save/load/resume run produce the same next iteration.

### Internal assignment representation

- Consider replacing pointer-valued `r_assign` and `q_assign` with stable `uint32_t` IDs only after model and
  evaluator work is stable.
- Do not renumber live IDs. Keep reusable object slots, a free-ID list, and dense per-locus active vectors.
- Retain linear swap-and-pop removal unless larger-cluster benchmarks show that position indexes repay their
  bookkeeping and memory cost.
- Consider `uint16_t` IDs only after measuring the maximum required live-ID capacity.

### Remaining performance work

- Profile maximization as removal, Viterbi scoring/backtracking, and reinsertion on fixed development and
  representative inputs.
- Record active R/Q counts, candidate counts, visited Q edges, cluster churn, allocation counts, peak memory,
  and cache misses.
- If numeric scoring dominates, benchmark cached transition terms, repeated emission terms, and integer
  `log(n)` values with precise invalidation rules.
- Run AddressSanitizer and UndefinedBehaviorSanitizer over hard, noisy, soft, sequential-init, PBWT-init,
  singleton-deletion, ID-reuse, and batched-maximization paths.
- Establish a final post-optimization baseline with exact commands, commit, compiler, threads, dimensions,
  output sizes, per-stage timings, and peak RSS.

### External imputation evaluation

- Compare DFCP with Beagle using the same reference, target, observed-marker mask, and windows.
- Compare pooled MAC-stratified r-squared and accuracy, runtime, and peak memory under matched conditions.
