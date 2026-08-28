# DFCP TODO

The immediate goal is to make maximization fast enough for routine development. Make and benchmark one
optimization at a time. Model persistence, resumable training, and cluster evaluation remain planned work,
but should wait until the fitting loop is fast enough to test reliably.

## Current state

- The chr20 pipeline windows aligned VCFs first, then bitpacks each generated window.
- `DFIP` probability output, `DFIE` imputation evaluation, and the imputation visualization are working.
- Reusing Viterbi message maps across sequences reduced maximization time by roughly half.
- The first window from a 5,000-window plan with overlap 32 is the development input: 4,904 reference
  haplotypes, 376 reference loci, 10 observed loci, and 366 masked loci. Hard PBWT training plus Viterbi
  imputation takes about seven seconds with four OpenMP threads after message reuse.
- The first window from a 500-window plan is the representative performance input: 4,904 reference
  haplotypes and 3,535 loci. A run took roughly four minutes before the current optimization baseline was
  fully recorded.

Do not compare timings from different inputs, modes, initialization methods, thread counts, compilers, or
commits.

## Decisions

### Training order

- Initialization constructs the first cluster graph using the prior parameter values.
- Each training iteration runs expectation, maximization, then ELBO evaluation.
- Do not add a separate expectation step between initialization and the loop. That extra update is not worth
  another full pass over the model.
- Sequential Viterbi initialization does not make the first maximization redundant: early sequences were
  inserted into a partial graph and may move after all reference sequences have been added.

### Assignments and cluster identity

- Keep both `r_assign` and `q_assign` during training. Removing a sequence requires its exact `R_l` and `Q_l`
  clusters; a `Q_l` cannot be recovered uniquely from its parent and child `R` clusters.
- Replace assignment pointers with stable `uint32_t` cluster IDs when live cluster storage is redesigned.
- Do not renumber live IDs. Use reusable object slots, a free-ID list, and dense per-locus active vectors.
- Delete from active vectors with swap-and-pop and store each cluster's active-vector position.
- Keep the hard-emission filtering supplied by `rs_by_emit`, but replace each unordered set with a dense
  vector during the storage redesign.
- Consider `uint16_t` IDs only after measuring the maximum live ID and cluster count per locus.

### Model persistence

- The first model format will use magic `DFCM` and contain only dimensions, emission mode, fitted parameters,
  `R` and `Q` clusters, and graph edges.
- Do not add a version number or reference/window identity while the format is still changing rapidly.
- Use file-local cluster IDs. The serialized representation must not depend on live in-memory IDs.
- Rebuild derived data such as active vectors and `rs_by_emit` after loading.
- A frozen inference model does not need assignments. A resumable checkpoint needs both `r_assign` and
  `q_assign`, enabled by an explicit file flag.
- Store fitted `Params` in `DFCM`; do not create a second required parameter file. Once model output exists,
  remove the large parameter arrays from stdout JSON.

## Active performance work

### 1. Lock down correctness and measurement

Status: next.

- Compare the reused-message implementation with the implementation before message reuse on one small,
  deterministic input in hard, noisy, and soft modes. Compare selected paths or final assignments, not only
  aggregate JSON.
- Confirm Viterbi imputation probabilities are unchanged.
- Run AddressSanitizer and UndefinedBehaviorSanitizer. Reused pointer-keyed maps retain stale keys, so code
  may dereference a key only when the current candidate traversal proves that cluster is live and its message
  was written during the current Viterbi call.
- Record initialization, expectation, maximization, ELBO, and imputation times for both fixed inputs.
- Profile maximization as removal, Viterbi scoring/backtracking, and reinsertion. Also record active `R` and
  `Q` counts, compatible hard-emission candidates, visited `Q` edges, cluster churn, retained map capacity,
  allocation counts, peak memory, and cache misses.

### 2. Remove redundant Viterbi work

Status: partially complete. Message maps are already reused across a maximization pass, sequential
initialization, and target imputation. Benchmark the remaining changes after phase 1.

Apply and benchmark these independently, in this order:

1. Remove existing-cluster `b_msgs`. An existing `Q_l` has one child `R_(l+1)`, so its message is the child's
   `a` message. Retain only the one new-`Q` message per transition.
2. Reuse the `2L-1` Viterbi path buffer across sequences.
3. In hard mode, skip a `Q` edge when its child `R` emission cannot match the next observed allele.
4. Benchmark `boost::unordered_flat_map` for the remaining `a` messages. There are currently no retained
   references or iterators whose invalidation should prevent its use.

Do not reserve maps inside `get_viterbi_clusters`. Reused maps retain capacity, and the earlier reserve test
was slower. Reconsider one-time reservation only if measured map growth identifies it as useful.

### 3. Replace hash-based live cluster storage

Status: pending after the smaller message changes.

- Add reusable stable-ID object slots and dense active `R` and `Q` vectors at each locus.
- Replace the pointer-keyed ownership map and convert both assignment arrays to `uint32_t` IDs.
- Replace `rs`, `qs`, and hard-emission unordered sets with dense vectors.
- Store a `Q` cluster's one parent and one child directly instead of using one-element vectors.
- Replace message maps with dense arrays indexed by stable ID or dense active position, whichever benchmarks
  better. Use an epoch or validity array where not every entry is overwritten, especially hard `a` messages.
- Preserve constant-time deletion, stable identity, and hard-emission filtering.
- Test sequential and PBWT initialization, all emission modes, validation splitting, singleton deletion,
  cluster-ID reuse, graph mutations, and temporary changes to `HP.N`.
- Measure assignment memory, total cluster memory, churn, and runtime before considering narrower IDs.

### 4. Optimize scoring only if profiling supports it

Status: pending.

- Separate time spent in container traversal and hashing from time spent in numeric scoring.
- If logarithms and `delta_Elogx` are material, precompute integer `log(n)` values and cache transition terms
  that remain fixed during one maximization pass.
- Invalidate only values affected by cluster size or degree changes during removal and reinsertion.
- Compute repeated emission terms once per locus or active cluster where possible.
- Avoid an `L x N` cache unless its measured speedup justifies its memory and cache cost.

### 5. Parallelize only independent or sufficiently large work

Status: pending after dense serial storage is measured.

- Parallelize frozen target imputation across target sequences first; it does not mutate the fitted graph.
- Keep training-sequence reassignment serial unless the inference algorithm is deliberately changed.
- If cluster scoring remains expensive, benchmark the loops in maximization and forward/backward using a
  persistent OpenMP region rather than starting a region at every locus.
- Preserve deterministic reductions and compare serial and parallel probability files byte for byte.

### 6. Establish the post-optimization baseline

Status: pending.

- Run hard, noisy, and soft training on both fixed inputs, followed by Viterbi and forward-backward
  imputation.
- Record exact commands, commit, compiler, thread count, dimensions, initialization, iteration count, output
  sizes, wall time, peak resident memory, and per-stage timings.
- Confirm stdout remains one small JSON document, diagnostics remain on stderr, probability files decode,
  evaluation succeeds, and reports rebuild from saved output.
- Set a concrete maximization-time target before beginning model persistence.

## Completion rule

Do not combine independent optimizations in one benchmark change. Retain a change only when it builds without
warnings, passes its focused correctness comparison, and improves the same representative input or provides
a clearly measured memory benefit worth its runtime cost.

## Later backlog

### Automated regression coverage

- Add a small synthetic VCF fixture that checks overlapping windows, `--n-generate`, the manifest, VCF record
  alignment, window-local `observed_loci.txt`, bitpacked alleles, and all output dimensions.
- Add focused `DFIP` and `DFIE` tests for exact decoding, fixed-point tolerance, atomic replacement,
  malformed headers, mismatched dimensions, truncation, trailing data, constant inputs, and the `0.5`
  accuracy threshold.
- Keep one deterministic end-to-end fixture covering hard, noisy, and soft modes and both imputation methods.
- Verify masked indexes are the complement of observed indexes and VCF `min(AC, AN-AC)` agrees with counts
  decoded from `ref.bin`.

### Frozen model save and load

- Specify the exact `DFCM` byte layout beside its reader and writer.
- Save fitted parameters, cluster records, and graph edges with file-local IDs; stream and atomically replace
  the output where practical.
- Load a frozen model, rebuild ownership and derived indexes, and validate all dimensions and graph links at
  the I/O boundary.
- Require imputation from an in-memory model and the same saved-then-loaded model to be byte-identical.
- Replace parameter arrays in stdout JSON with the model output once this path is stable.

### Resumable training checkpoints

- Extend `DFCM` with an explicit assignments-present flag only after frozen loading works.
- Store and restore both assignment arrays and every state value required to continue training exactly.
- Verify that an uninterrupted run and a save/load/resume run produce the same next iteration.

### Cluster and tree evaluation

- Restore the inactive evaluator against the saved-model representation and current tree APIs.
- Test exact clades, sibling-split clusters, singleton and full-sample clusters, dominated ancestors, random
  small trees against brute force, and validation runs with held-out leaves.
- Restore excess parsimony, importance-weighted cluster-to-clade IoU, adjacent-locus IoU, emission IoU,
  purity, and tree visualization only after the focused tests pass.

### External imputation evaluation

- Compare DFCP with Beagle on the same reference, target, observed-marker mask, and windows.
- Compare per-locus and MAC-stratified r-squared and accuracy, and record runtime and peak memory under
  matched conditions.
