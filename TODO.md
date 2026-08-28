# DFCP TODO

The immediate goal is to make maximization fast enough for routine development. Make and benchmark one
optimization at a time. Model persistence, resumable training, and cluster evaluation remain planned work,
but should wait until the fitting loop is fast enough to test reliably.

## Current state

- The chr20 pipeline windows aligned VCFs first, then bitpacks each generated window.
- `DFIP` probability output, `DFIE` imputation evaluation, and the imputation visualization are working.
- Reusing `ViterbiBuffers` and indexing inference messages by reusable cluster IDs removed allocation and
  hashing from Viterbi and forward-backward scoring.
- The first window from a 5,000-window plan with overlap 32 is the development input: 4,904 reference
  haplotypes, 376 reference loci, 10 observed loci, and 366 masked loci. Hard PBWT training plus Viterbi
  imputation takes about 4.2 seconds with four OpenMP threads after the current Viterbi optimizations.
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
- Active vectors use swap-and-pop after a linear pointer search. At about 7.3 R clusters per locus, storing
  positions made maximization 5.6% slower and added about 1.3 MB peak RSS. On a fixed 376-locus input,
  linear and indexed removal remained effectively tied at PBWT match lengths 10 and 20, with about 13.9 and
  31.9 final R clusters per locus; indexed removal still used more memory, so retain the simpler search.
- Keep the hard-emission filtering supplied by the dense `rs_by_emit` vectors.
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
- Run AddressSanitizer and UndefinedBehaviorSanitizer. Reused dense message entries may be read only when
  current-candidate traversal proves that the entry was written during the current inference pass.
- Record initialization, expectation, maximization, ELBO, and imputation times for both fixed inputs.
- Profile maximization as removal, Viterbi scoring/backtracking, and reinsertion. Also record active `R` and
  `Q` counts, compatible hard-emission candidates, visited `Q` edges, cluster churn, dense message capacity,
  allocation counts, peak memory, and cache misses.

### 2. Remove redundant Viterbi work

Status: complete.

Apply and benchmark these independently, in this order:

1. Reuse the `2L-1` Viterbi path buffer across sequences. Complete: this removed about 54,000 allocation and
   free pairs from the development run, but changed wall and maximization time by less than 0.2%.
2. Remove existing-cluster `b_msgs`. Complete: existing `Q_l` clusters now read their sole child's `a`
   message directly, leaving one `new_b_msg` per transition. Development-window maximization became 22.6%
   faster and wall time became 22.2% faster.
3. Skip hard-mode `Q` edges whose child emission cannot match the next observed allele. Complete: the early
   skip and direct current-message lookup made development-window maximization 4.4% faster and wall time
   4.3% faster.
4. Benchmark `boost::unordered_flat_map` for the remaining `a` messages. Complete, not retained: ten
   alternating development-window runs improved mean maximization time by 1.2% and median time by 0.6%,
   with no memory reduction. This is too small to carry into the cluster-ID redesign.

### 3. Replace hash-based live cluster storage

Status: in progress. Clusters now have reusable `uint32_t` IDs, and inference messages use dense ID-indexed
vectors. Five alternating development-window runs reduced mean maximization time by 16.4%. On the
representative window, one run reduced maximization from 91.6 to 63.7 seconds; excluding one baseline
outlier gives a 26.6% per-iteration reduction. Dense forward-backward messages reduced mean representative
imputation time by 33.2%. Replacing active hash sets with dense vectors and linear removal reduced mean
representative maximization time another 20.6%, from 5,628 to 4,467 ms per iteration. Hard, noisy, and soft
probability files were byte-identical across the applicable changes. ID-indexed ownership and direct Q links
then reduced default representative maximization from 4,467 to 4,139 ms per iteration and wall time from
about 59.1 to 54.8 seconds, while lowering peak RSS by about 3.7 MB.

- The pointer-keyed ownership map has been replaced by reusable ID-indexed `unique_ptr` slots. At PBWT match
  length 20 this reduced mean maximization and wall time by 2.4% and peak RSS by about 1.1 MB.
- Convert both assignment arrays to `uint32_t` IDs.
- `rs`, `qs`, and the hard-emission index are dense vectors. Linear removal beat stored vector positions on
  the representative window, and every training iteration now records the mean R-cluster count per locus.
- Q clusters store their one parent and one child directly. At PBWT match length 20 this reduced mean
  maximization and wall time by another 6.8%, with effectively unchanged peak RSS.
- Viterbi and forward-backward message maps have been replaced by dense arrays indexed by stable ID.
  Current-candidate traversal guarantees that every accessed hard-mode message was overwritten during the
  current inference pass.
- Preserve stable identity and hard-emission filtering.
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
