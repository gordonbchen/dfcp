# DFCP TODO

## Current state

- `scripts/fsc_sim/run.sh` runs fastsimcoal with an explicit seed and prepares its output.
- `scripts/fsc_sim/prep_data.py` converts one haploid `.gen` file into aligned `ref.bin` and
  `variant_pos.txt` files.
- `impute` can train with only a reference file. The target, observed-loci, and probability files are one
  optional positional group.
- `--output_r_assign FILE` writes final reference R assignments in the `DFRA` format without requiring
  target imputation.
- `eval_clusters` reads reference alleles, `DFRA` R assignments, variant positions, and fastsimcoal true
  trees. It reports partition continuity, cluster count, purity, excess parsimony, and cluster-to-clade IoU.
- The complete 100-haplotype, 13,624-locus fastsimcoal fixture trains and evaluates successfully.

## Next cluster-evaluation work

- Review which scalar metrics should remain in the stable `eval_clusters` JSON interface.
- Update or replace the legacy `scripts/parsimony.py`, which still targets the removed text-sequence
  executable interface.
- Add a small permanent fixture for sparse cluster IDs, exact clades, split clades, singleton clusters,
  full-sample clusters, repeated positions, and recombination-boundary selection.
- Compare random small-tree results with brute-force descendant sets and parsimony enumeration.
- Decide whether tree visualization belongs in `eval_clusters` or a separate visualization tool.

## Later backlog

### Frozen models and separate training

- Reconsider splitting fitting and inference only when repeated imputation from one fitted model makes the
  additional format and executable worthwhile.
- Define `DFCM` beside its reader and writer. Store dimensions, emission mode, inference parameters, and
  file-local R/Q records. Each Q record can identify its parent and child R; rebuild R adjacency and derived
  indexes while loading.
- Require in-memory and saved-then-loaded models to produce byte-identical probabilities in hard, noisy,
  and soft modes before changing the command-line interface.
- If the split is retained, add `train REF MODEL [options]` and change `impute` to load `MODEL`.

### Resumable training checkpoints

- Extend model serialization with assignments only after frozen inference works.
- Store both R and Q assignments and all mutable state needed to resume training exactly.
- Verify that uninterrupted and save/load/resume runs produce the same next iteration.

### Maintained script interfaces

- Update `scripts/dfcp.py`, `scripts/init.py`, and `scripts/tune.py` to the current binary-sequence interface
  where they remain useful.
- Keep generated reports rebuildable without rerunning DFCP.

### Automated imputation regression coverage

- Add a small VCF fixture covering overlapping windows, `--n-generate`, manifest rows, aligned records,
  observed-locus indexes, bitpacked alleles, and output dimensions.
- Add focused `DFIP` and pooled-evaluator tests for fixed-point decoding, overlap ownership, minor-allele
  orientation, mismatched dimensions, constant bins, and the `0.5` hard-call threshold.
- Keep one deterministic end-to-end fixture covering hard, noisy, and soft modes and both imputation methods.
- Verify masked indexes are the complement of observed indexes and VCF `min(AC, AN-AC)` agrees with decoded
  reference counts.

### Internal assignment representation

- Consider replacing pointer-valued `r_assign` and `q_assign` with stable `uint32_t` IDs only after the
  evaluator is stable.
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
