# Target-specific reference filtering

## First experiment: whole-window weighted top K

For each target haplotype, rank the reference haplotypes using only target
alleles at observed loci. Keep the top `K_ref` references, then construct PBWT
blocks, fit DFCP, and impute that target using the filtered panel. Do not use
masked target truth in selection.

At observed locus `l`, let `t_l` be the target allele (0 or 1), `r_l` the
reference allele, and `p_l(a)` the frequency of allele `a` in the **full,
unfiltered reference panel**. Define

```text
score(r, t) = sum over observed l of 1[r_l = t_l] * (-log p_l(t_l))
```

Thus a matching 0 or 1 contributes `-log` of that allele's frequency;
a mismatch contributes **zero**. `p_l(t_l)` is not always the minor
allele frequency: use the frequency of the allele actually carried by the
target. A rare matching allele is more informative than a common matching
allele. The score is a ranking heuristic, not a probability that the entire
reference haplotype was copied.

Estimate frequencies as `allele_count / N_ref` from all reference haplotypes
before filtering. A reference cannot match an allele absent from the panel,
so the score never needs `log(0)`; a monomorphic match has weight zero. This
first experiment has no frequency smoothing or weight cap. Break equal-score
ties by original reference index for reproducibility. If `K_ref >= N_ref`,
keep the whole panel. If the top K omit every carrier of an observed target
allele, append the highest-scoring reference that carries it. This is needed
for block compatibility and makes K a lower bound on panel size rather than
an exact size. Record the actual selected size for every target.

## Why top K rather than a fixed score threshold?

The 64-window observed-allele histogram is in `output/ref_match_hist.html`;
`scripts/ref_match_viz.py` regenerates it. A fixed percentage of matching
sites gives very different panel sizes across targets. For example, 80%
agreement retained a median of 1,316 out of 4,904 references but retained
none for one target-window pair; 85% retained a median of 168 and none for
182 target-window pairs. Top K gives a predictable approximate fit size and avoids an
empty panel. The histogram uses unweighted agreement, so it does not itself
choose `K_ref` for the weighted score.

## 64-window experiment

The experiment used `data/1000g_phase3_v5b/windows_64`, emission
initialization, `--block_max_k 1024`, `--max_batch_size 64`, one training step,
eight window workers, and `OMP_NUM_THREADS=1`. The unfiltered comparison is
`individual_k1024_step1`, which used the same fitting settings. Every run
evaluated the same 64 windows and the same 1,665,200 masked loci. Wall times
are measured for the imputation runner, excluding `eval_impute`.

| Reference selection | Wall time | Speedup vs full | Mean blocks/target | r² MAC 1 | r² MAC 10 | r² MAC 100 | r² MAC 500 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Full 4,904 | 12m56s | 1.00× | 36.38 | 0.0901 | 0.5227 | 0.8643 | 0.9650 |
| Top 256 + coverage | 54s | 14.33× | 5.23 | 0.0727 | 0.5453 | 0.8395 | 0.9377 |
| Top 512 + coverage | 1m48s | 7.22× | 4.94 | 0.0781 | 0.5680 | 0.8498 | 0.9445 |
| Top 1,024 + coverage | 3m18s | 3.93× | 4.77 | 0.0815 | 0.5719 | 0.8526 | 0.9368 |
| Top 2,048 + coverage | 6m04s | 2.13× | 18.26 | 0.1113 | 0.6019 | 0.8848 | 0.9570 |
| Top 3,000 + coverage | 8m16s | 1.57× | 25.53 | 0.1010 | 0.5644 | 0.8812 | 0.9550 |

The coverage rule added 74, 30, 11, 2, and 0 references in total across
6,656 target fits at the five filtered Ks, respectively. Top 2,048 is the
best tested speed–accuracy tradeoff: it improves low-MAC r² over the full
panel at these representative bins, with a slight loss at MAC 500. Top 3,000
is slower and less accurate than top 2,048 at every representative bin above.
Across broader MAC ranges, top 3,000 has a small average advantage above
MAC 100 but a larger deficit below MAC 50; those averages are weighted means
of per-MAC r² values, not pooled r².
Filtering also changes the block grid: at `K_ref <= block_max_k`, the pattern
cap cannot trigger, so target compatibility largely determines block lengths.
The measured effect is therefore not solely from fewer training haplotypes.
All r² values are pooled within each exact MAC bin; selected bins above are
not an aggregate r² over a broader frequency range.

The curves are in `output/ref_filter_compare.html`. Each run has per-window
probabilities, JSON diagnostics, and logs under
`data/1000g_phase3_v5b/windows_64/window_*/impute/ref{K}_b1024_step1.*`;
evaluations are in `output/ref{K}_b1024_step1.tsv`.

To regenerate one run and its evaluation:

```bash
OMP_NUM_THREADS=1 scripts/1000g_phase3_v5b/impute.sh \
  data/1000g_phase3_v5b/windows_64 \
  --name ref2048_b1024_step1 --n-parallel 8 --init emission \
  --block_max_k 1024 --ref_top_k 2048 \
  --max_batch_size 64 --max_train_steps 1
./build/eval_impute data/1000g_phase3_v5b/windows_64 \
  ref2048_b1024_step1 output/ref2048_b1024_step1.tsv
```

## Verification

- `scripts/test_ref_filter.cpp` checks matching 0s and 1s, tie order,
  a monomorphic reference locus, observed-allele coverage, and `K_ref >= N_ref`.
- A small end-to-end fixture produced byte-identical probability files with
  filtering disabled and with `K_ref = N_ref`.
- Each of the five filtered runs completed 64 windows and 6,656 target fits.
  All probability-file dimensions and JSON target counts were checked. The
  six evaluation TSVs have identical MAC bins and prediction counts.

## Risks and later alternatives

Filtering changes the data used to form PBWT blocks, cluster assignments,
emission frequencies, and fitted parameters. It can therefore alter cluster
structure, not just runtime. A high whole-window score can also hide a poor
local match, while a reference useful over one tract can rank poorly over the
whole window. Check whether selected references still provide compatible
patterns at observed target sites; do not silently accept failed construction.

If whole-window top K hurts accuracy, the next experiment is to select some
references from each part of the window, or to use PBWT-based local tract
matches. Keep those as separate experiments rather than adding them to the
first implementation.
