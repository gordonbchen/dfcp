# DFCP TODO

## Active: greedy PBWT-block DFCP

### Model and representation

- [x] Add opt-in `--block_max_k K`. Starting at each boundary, extend through the longest SNP prefix having
  at most `K` exact reference patterns. With target inputs, also require every target haplotype to retain at
  least one reference sequence matching all of its observed alleles in the candidate block. Cut immediately
  before the first SNP that violates either monotone constraint. A zero maximum retains the SNP-locus model.
- [x] Maintain target-compatible reference sets as packed bitsets. Reset them at each cut and never inspect
  masked target truth. A single observed target allele absent from the reference remains an input error; the
  reference-polymorphism preparation step prevents this on the maintained 1000 Genomes panel.
- [x] Select boundaries in one `O(N)`-working-memory streaming PBWT pass, then run a second streaming pass to
  write dense sequence-major labels. This avoids full `N * L` PBWT tables and avoids retaining block-major
  labels merely to transpose them.
- [x] Store source-SNP offsets and reference counts for every emission. Report mean block length, mean
  emissions per block, the singleton block fraction, and counts of maximum-K and compatibility cuts.
- [x] Store model observations as dense block-label integers in sequence-major order. Retain one reference per
  label so SNP alleles can be recovered without duplicating every pattern. Use compact offsets for the ragged
  `K_b` alphabet.
- [x] Keep source-SNP and model-block coordinates distinct. `HP.L`, parameters, assignments, and messages use
  model blocks; target masks, probability output, and evaluation remain in source-SNP space.
- [x] Initialize block R clusters from block labels and Q clusters from distinct adjacent R pairs. The PBWT
  radius initializer remains a SNP-mode heuristic and is not applied to categorical block labels.
- [x] Check random small references and target mosaics against brute-force pattern equality, `K_b`,
  compatibility, emission counts, and boundary maximality. Recheck exact SNP-mode regressions.

### Target inference and outputs

- [x] Convert each target's observed SNPs into a compatible-label set for every block: label `e` is compatible
  exactly when its representative reference haplotype agrees at every observed SNP in that block. Fail with a
  useful error if a target/block has no compatible reference label; polymorphism at each SNP alone does not
  mathematically guarantee that a joint multi-marker pattern exists.
- [x] Generalize exact-emission Viterbi and forward-backward candidates from one label (or all labels)
  to the compatible-label set. Preserve the proposed-new-R state: its emission term is the log-sum of the
  categorical predictive masses of compatible labels, not an unweighted compatibility indicator.
- [x] Emit SNP-level allele probabilities. An existing R state contributes the deterministic allele of its
  block label. A proposed-new-R state contributes a mixture over compatible labels weighted by the same
  categorical predictive probabilities used in its emission message. Observed target SNPs remain constraints;
  only the complement of `observed_loci.txt` is written, in increasing original SNP order.
- [x] Expand each block R assignment across its SNPs when writing `DFRA`, preserving the current `N x L_snps`
  evaluator contract. IDs remain equality labels; repeating one ID throughout a block deliberately makes
  boundaries visible to the existing tract metrics without changing `eval_clusters`.
- [x] Add JSON diagnostics needed to interpret compression and cost: maximum K, SNP and block counts, mean
  block length, mean emissions per block, occupancy and singleton summaries, and cut reasons.

### Correctness checks

- [x] Compare random greedy boundaries and label equivalence classes with brute-force reference patterns.
- [x] Verify the zero-maximum SNP path retains byte-identical fitted assignments and imputation probabilities.
- [x] Cover target blocks with zero, one, and several observed SNPs; multiple compatible labels; incompatible
  joint observations; a masked SNP in the same block as an observed SNP; and a short last block.
- [x] Build with `./build.sh`, run sanitizer builds over the new coordinate conversions, and confirm stdout is
  one JSON object while diagnostics remain on stderr.

### Benchmarks and decision gates

- [ ] Benchmark caching target-compatible block emissions and R candidates across identical observed target
  patterns during imputation. `BlockObs` currently rebuilds them for every target sequence; keep a cache only
  if its memory cost is justified by an end-to-end speedup.

- [x] Record a fresh SNP-mode baseline before judging the new path. Fix commit, compiler, `OMP_NUM_THREADS`,
  batch size, train steps, initialization settings, wall time, stage timings, peak RSS, and output sizes.
- [x] On the fastsimcoal fixture, compare maximum K values `4, 8, 16, 32, 64, 128` with the same priors and
  training budget. Expand assignments to SNPs and report `mean_adj_iou`, `mean_clusters`, excess parsimony,
  clade IoU, exact tract loci/bp, fit time, and peak RSS. Include initialization-only and trained runs so the
  model update is distinguishable from the PBWT partition itself.
- [x] Across all 1000 Genomes windows, compare maximum K values `4, 8, 16, 32, 64, 128, 256` using pooled
  MAC-stratified imputation r-squared, initialization/training/imputation time, and peak RSS. Use
  the same target, mask, overlap ownership, priors, batch size, threads, and train steps for every run.
- [x] Scale the 1000 Genomes reference by joining prepared windows or generating larger windows only after
  single-window correctness. Test increasing SNP counts until memory or time is clearly limiting; do not
  start with a full chromosome. Report throughput per reference SNP and per model block as well as peak RSS.
- [x] Treat maximum K as a quality/performance frontier, not a compression-only win. A useful default
  must reduce model loci and wall time materially without a large loss in rare-variant r-squared or clade IoU.
  Keep the SNP path if no maximum meets both sides of that gate.

Baseline recorded 2026-09-16 from commit `6861140`, GCC 16.2.1, `OMP_NUM_THREADS=1`, exact emissions,
PBWT initialization, and GNU time 1.10:

- fastsimcoal (`N=100`, `L=13,624`, match length 200): initialization-only was 0.49 s wall and 505 MiB peak
  RSS. Three ME steps were 22.81 s wall and 533 MiB peak RSS; maximization totaled 20.54 s. The expanded
  SNP-locus evaluation gave mean adjacent IoU 0.9861, 29.39 mean clusters, 2.238 excess parsimony, clade IoU
  0.8698, and mean exact-tract span 29,660 bp.
- 1000 Genomes window 0007 (`N_ref=4,904`, `L=3,535`, `N_target=104`, 216 observed target SNPs, match length
  20, batch size 4): initialization-only fitting plus forward-backward imputation was 5.08 s wall and 493 MiB
  peak RSS (`t_init=502 ms`, `t_impute=4,544 ms`). One ME step was 106.46 s wall and 493 MiB peak RSS;
  `t_max=96,859 ms` and `t_impute=8,856 ms`.
- Existing initialization-only `pbwt20` probabilities across all 500 prepared 1000 Genomes windows evaluate in
  9.78 s. These predate this work and are a quality reference, not a fresh fit: MAC-1 r-squared is
  `2.13e-6`, MAC-10 is `0.00587`, and MAC-20 is `0.0351`.

### Greedy-block benchmark results

An earlier equal-width prototype was removed: block sizes `2,4,8,16,64` rejected respectively
`1,3,5,21,118` of 500 windows because per-SNP polymorphism does not ensure a matching joint target pattern.
The target-aware greedy construction produced all 500 windows at every tested maximum. Runs used one ME step,
one OpenMP thread per process, and 12 simultaneous processes.
`singleton assign` is the fraction of all reference haplotype/block assignments whose emission contains only
that haplotype; it is more representative of exposure than the unweighted singleton-emission fraction.

| max `K_b` | mean blocks/window | mean SNPs/block | mean haps/emission | singleton assign | compatibility cuts | wall s | max RSS MiB | MAC-1 r2 | MAC-10 r2 | MAC-20 r2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 1179.7 | 3.00 | 1232.4 | 0.027% | 3 | 54.55--58.56 | 123.2 | 0.0000393 | 0.000643 | 0.00308 |
| 8 | 512.9 | 6.89 | 617.8 | 0.064% | 14 | 33.17--41.14 | 60.6 | 0.0000540 | 0.00174 | 0.00757 |
| 16 | 245.0 | 14.43 | 309.5 | 0.141% | 17 | 25.68--31.43 | 35.9 | 0.000126 | 0.00438 | 0.0188 |
| 32 | 122.9 | 28.77 | 155.2 | 0.301% | 50 | 25.11 | 24.4 | 0.000873 | 0.0100 | 0.0562 |
| 64 | 63.4 | 55.72 | 77.9 | 0.634% | 118 | 28.61 | 23.7 | 0.00648 | 0.0500 | 0.156 |
| 128 | 33.5 | 105.66 | 39.5 | 1.330% | 261 | 36.14 | 23.8 | 0.0190 | 0.162 | 0.340 |
| 256 | 18.0 | 196.07 | 20.4 | 2.757% | 571 | 46.00 | 23.7 | 0.0406 | 0.290 | 0.476 |

Across the panel, initialization/maximization/imputation totals in seconds were respectively
`67.7/357.0/133.2`, `31.7/172.4/121.5`, `15.5/89.5/139.6`, `8.2/46.5/184.7`, `4.5/24.5/256.4`,
`2.6/14.0/359.6`, and `1.8/8.4/484.0` as maximum K increased. The fastest end-to-end region is therefore
`K_b=16--32`; larger maxima trade more forward-backward work and reference specificity for much stronger
imputation r-squared. A repeat of the
three smallest maxima after cleanup produced the wall-time ranges shown above; maximum RSS was stable, while
concurrent stage timings varied by 8--25%, so small timing differences should not be overinterpreted.

On fastsimcoal, three-step results were:

| max `K_b` | blocks | mean SNPs/block | singleton assign | wall s | RSS MiB | adjacent IoU | mean clusters | excess parsimony | clade IoU | tract bp |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 3955 | 3.44 | 0.66% | 0.73 | 23.4 | 0.8563 | 4.23 | 0.185 | 0.9057 | 2012 |
| 8 | 1604 | 8.49 | 1.81% | 0.23 | 17.4 | 0.9286 | 8.58 | 0.507 | 0.8605 | 6547 |
| 16 | 691 | 19.72 | 4.90% | 0.11 | 13.6 | 0.9663 | 17.70 | 0.745 | 0.8341 | 18866 |
| 32 | 274 | 49.72 | 13.64% | 0.07 | 11.2 | 0.9871 | 36.50 | 0.397 | 0.8443 | 65443 |
| 64 | 67 | 203.34 | 43.01% | 0.02 | 7.9 | 0.9975 | 72.46 | 0.025 | 0.9211 | 646869 |
| 128 | 1 | 13624 | 98.00% | <0.01 | 5.4 | 1.0000 | 99.00 | 0 | 1.0000 | 9997970 |

The `K_b=128` result confirms the overfitting concern directly: the entire simulated chromosome becomes one
block with 99 patterns for 100 reference haplotypes, and 98% of assignments are singletons. A raw `K_b` cap is
therefore a real regularizer, not just a memory limit.

A target-aware no-window chromosome-20 initialization and forward-backward run at `K_b=64` used 30,604 blocks
(mean 56.1, range 3--85 SNPs), 1,939,659 emissions, and 122 compatibility cuts. It completed all 104 target
haplotypes in 287.9 s wall (`t_impute=259.6 s`) at 5.15 GiB peak RSS and wrote 331 MiB. Full-chromosome
r-squared at MAC `1,10,20` was `0.00186,0.0555,0.164`. Thus adaptive construction makes unwindowed imputation
valid, but windowed parallel execution is substantially faster. Against windowed one-step `K_b=64`, the
no-window initialization was lower at MAC 1 (`0.00186` versus `0.00648`) and similar at MAC 10/20.

### PBWT comparison

Full-panel PBWT initialization-only imputation at radii `10`, `20`, and `50` gave pooled r-squared at MAC
`1/10/20` of `1.71e-6/0.00353/0.0158`, `2.13e-6/0.00587/0.0351`, and
`5.49e-6/0.00816/0.0743`, respectively. Larger PBWT radii improve common-variant imputation but remain below
the corresponding larger greedy blocks in the table above.

Current single-thread measurements on 1000 Genomes window 0007 (`N_ref=4,904`, 3,535 SNPs, 104 targets) used
the same initialization-only PBWT setup and one ME step for greedy blocks:

| method | parameter | model loci | wall s | peak RSS MiB | init ms | ME ms | impute ms |
|---|---:|---:|---:|---:|---:|---:|---:|
| PBWT | radius 10 | 3,535 | 1.71 | 430 | 389 | 0 | 1309 |
| PBWT | radius 20 | 3,535 | 4.22 | 468 | 397 | 0 | 3806 |
| PBWT | radius 50 | 3,535 | 13.64 | 636 | 510 | 0 | 13061 |
| greedy | `K_b=4` | 1,202 | 0.78 | 122 | 91 | 506 | 118 |
| greedy | `K_b=16` | 268 | 0.36 | 35 | 17 | 124 | 171 |
| greedy | `K_b=64` | 80 | 0.41 | 18 | 7 | 49 | 310 |
| greedy | `K_b=256` | 29 | 0.62 | 14 | 4 | 31 | 539 |

Current fastsimcoal three-step runs (`N=100`, 13,624 loci, one thread) show that PBWT's reference-radius
initialization is much more memory-intensive than greedy blocks, while its cluster quality lies between the
small and large greedy maxima:

| method | parameter | wall s | RSS MiB | adjacent IoU | mean clusters | excess | clade IoU | tract bp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| PBWT | radius 50, init | 0.24 | 237 | 0.9797 | 47.97 | 0.162 | 0.8717 | 48,533 |
| PBWT | radius 100, init | 0.26 | 305 | 0.9899 | 64.35 | 0.047 | 0.9037 | 140,534 |
| PBWT | radius 200, init | 0.31 | 375 | 0.9952 | 77.75 | 0.010 | 0.9326 | 449,504 |
| PBWT | radius 200, 3 steps | 18.16 | 399 | 0.9861 | 29.39 | 2.238 | 0.8698 | 29,660 |
| greedy | `K_b=4`, 3 steps | 0.55 | 22 | 0.8563 | 4.23 | 0.185 | 0.9057 | 2,012 |
| greedy | `K_b=16`, 3 steps | 0.09 | 13 | 0.9663 | 17.70 | 0.745 | 0.8341 | 18,866 |
| greedy | `K_b=64`, 3 steps | 0.02 | 7 | 0.9975 | 72.46 | 0.025 | 0.9211 | 646,869 |
| greedy | `K_b=128`, 3 steps | <0.01 | 5 | 1.0000 | 99.00 | 0 | 1.0000 | 9,997,970 |

Generated comparisons are `output/benchmarks/imputation/comparison.html` for all 500-window imputation
curves and `output/fsc_blocks/comparison.html` for clade-time and tract-length densities. The latter includes
PBWT radii `50/100/200`, PBWT radius 200 after three steps, and greedy `K_b=4...128` after three steps.

## Backlog

### Cluster evaluation

- Review which scalar metrics should remain in the stable `eval_clusters` JSON interface.
- Add a small permanent fixture for sparse cluster IDs, exact clades, split clades, singleton clusters,
  full-sample clusters, repeated positions, and recombination-boundary selection.
- Compare random small-tree results with brute-force descendant sets and parsimony enumeration.
- Consider a lagged pairwise-IoU decay curve if exact tract lengths and adjacent-locus IoU do not adequately
  describe partition persistence.
- Keep the Beagle 4 DAG-edge baseline separate from Beagle 5: its local DAG partitions all haplotypes, whereas
  Beagle 5 composite-reference states do not.

### Frozen models and separate training

- Reconsider splitting fitting and inference only when repeated imputation from one fitted model makes the
  additional format and executable worthwhile.
- Greedy compatibility-aware boundaries currently use the supplied panel's observed alleles, so the fitted
  model is intentionally target-panel-specific. A separately trained reusable reference model would need to
  store its block offsets and can still encounter an incompatible future target. Before splitting training and
  imputation, choose explicitly between boundaries learned from a representative/calibration target panel, a
  conservative shared-marker-mask rule, or a later nonzero target observation model. Do not silently weaken
  the current exact-compatibility guarantee.
- Define `DFCM` beside its reader and writer. Store dimensions, inference parameters, and
  file-local R/Q records. Each Q record can identify its parent and child R; rebuild R adjacency and derived
  indexes while loading.
- Require in-memory and saved-then-loaded models to produce byte-identical probabilities before changing the
  command-line interface.
- If the split is retained, add `train REF MODEL [options]` and change `impute` to load `MODEL`.

### Resumable training checkpoints

- Extend model serialization with assignments only after frozen inference works.
- Store both R and Q assignments and all mutable state needed to resume training exactly.
- Verify that uninterrupted and save/load/resume runs produce the same next iteration.

### Automated imputation regression coverage

- Add a small VCF fixture covering overlapping windows, `--n-generate`, manifest rows, aligned records,
  observed-locus indexes, bitpacked alleles, and output dimensions.
- Add focused `DFIP` and pooled-evaluator tests for fixed-point decoding, overlap ownership, minor-allele
  orientation, mismatched dimensions, constant bins, and the `0.5` hard-call threshold.
- Keep one deterministic end-to-end fixture covering both imputation methods.
- Verify masked indexes are the complement of observed indexes and VCF `min(AC, AN-AC)` agrees with decoded
  reference counts.

### Internal assignment representation

- Consider replacing pointer-valued `r_assign` and `q_assign` with stable `uint32_t` IDs only after the
  evaluator
  is stable.
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
- Run AddressSanitizer and UndefinedBehaviorSanitizer over sequential-init, PBWT-init,
  singleton-deletion, ID-reuse, and batched-maximization paths.
- Establish a final post-optimization baseline with exact commands, commit, compiler, threads, dimensions,
  output sizes, per-stage timings, and peak RSS.

### External imputation evaluation

- Compare DFCP with Beagle using the same reference, target, observed-marker mask, and windows.
- Compare pooled MAC-stratified r-squared, runtime, and peak memory under matched conditions.
