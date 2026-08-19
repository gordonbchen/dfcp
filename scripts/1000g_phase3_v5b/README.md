# 1000 Genomes Phase 3 v5b preprocessing

Run the complete chromosome 20 pipeline from the repository root:

```bash
THREADS=8 scripts/1000g_phase3_v5b/run.sh
```

Downloads are skipped when the destination already exists. To download every
raw file again, use:

```bash
THREADS=8 scripts/1000g_phase3_v5b/run.sh --force
```

The five stages and their outputs are:

1. `download.sh` downloads the Phase 3 v5b VCF, its index, the population
   panel, the Illumina Omni2.5 BED, and the HapMap Phase II GRCh37 genetic-map
   archive into `data/1000g_phase3_v5b/`. It extracts the PLINK chr20 map for
   the window report.
2. `biallelic_snvs.sh` retains PASS biallelic SNVs and writes
   `biallelic_snvs.vcf.gz`, its CSI index, and
   `biallelic_snvs.svtype_counts.tsv` into `biallelic_snvs/`.
3. `split_ref_target.py` selects two target individuals per population and
   writes the reference and target VCFs, indexes, and sample lists into
   `ref_target/`.
4. `mask_target.sh` writes the chromosome 20 Omni BED, observed target VCF,
   withheld target truth VCF, and indexes into `mask_target/`.
5. `dfcp_prep.sh` streams each VCF directly into locus-major, 64-bit-padded
   DFCP inputs `ref.bin`, `target_observed.bin`, and
   `target_masked_true.bin` in `dfcp_prep/`. It also writes
   `observed_loci.txt`, which maps observed target columns to reference loci.
   The C++ reader performs the packed transpose.

After preparing the data, write only the first small overlapping window with:

```bash
python3 scripts/1000g_phase3_v5b/window.py \
  --output-dir data/1000g_phase3_v5b/windows_18000_o32 \
  --n-windows 18000 --overlap 32 --first-only
```

This produces a 128-locus first window for the current chr20 inputs. Omit
`--first-only` to produce all 18,000 windows. The script copies locus-major
packed records without decoding their alleles and writes `windows.tsv` with
the chromosome-wide bounds and counts.

Reproduce the interactive window-selection report with:

```bash
.venv/bin/python scripts/1000g_phase3_v5b/window_viz.py --output docs/window.html
```

Use `window_viz.py --help` for input overrides. `--physical-only` builds the
same report in Mb without a genetic map.

Train, impute, evaluate, and visualize the first window with:

```bash
window_dir=data/1000g_phase3_v5b/windows_18000_o32/window_0000
./build/dfcp "$window_dir/ref.bin" "$window_dir/target_observed.bin" \
  "$window_dir/observed_loci.txt" "$window_dir/soft_fwd_bkwd.probs.bin" \
  --mode soft --init pbwt --viterbi_impute 0 > "$window_dir/soft_fwd_bkwd.json"
./build/eval_impute "$window_dir/soft_fwd_bkwd.probs.bin" \
  "$window_dir/target_masked_true.bin" "$window_dir/soft_fwd_bkwd.eval.bin"
.venv/bin/python scripts/impute_viz.py "$window_dir/soft_fwd_bkwd.eval.bin" \
  "$window_dir/ref.bin" "$window_dir/observed_loci.txt" \
  --output "$window_dir/soft_fwd_bkwd.html"
```

## Recorded first-window baseline

Measured on 2026-08-18, window 0 had 4,904 reference haplotypes, 104 target
haplotypes, 128 reference loci, 4 observed loci, and 124 masked loci. Creating
it took 0.048 seconds. Every probability file was 25,804 bytes and every
evaluation file was 1,000 bytes, exactly matching their declared dimensions.

| Mode | Imputation | ME iterations | Total time | Imputation | Mean r² | Defined r² | Accuracy |
|---|---|---:|---:|---:|---:|---:|---:|
| hard | Viterbi | 10 | 3.696 s | 9 ms | -1 | 0 / 124 | 0.991160 |
| hard | forward-backward | 10 | 3.798 s | 27 ms | 0.041987 | 6 / 124 | 0.991160 |
| soft | Viterbi | 8 | 2.146 s | 5 ms | 0.028945 | 2 / 124 | 0.991160 |
| soft | forward-backward | 8 | 2.145 s | 14 ms | 0.090874 | 8 / 124 | 0.991160 |

A repeated soft forward-backward run used 25,464 KiB peak resident memory and
2.157 seconds wall time. The low number of defined r² values is expected for
this exceptionally small, sparse first window: many truth or predicted vectors
are constant across the 104 targets. Synthetic tests separately verify defined
and undefined r² behavior.

Run stages individually with their defaults if needed:

```bash
scripts/1000g_phase3_v5b/download.sh
scripts/1000g_phase3_v5b/biallelic_snvs.sh
python3 scripts/1000g_phase3_v5b/split_ref_target.py --seed 0 --threads 8
THREADS=8 scripts/1000g_phase3_v5b/mask_target.sh
scripts/1000g_phase3_v5b/dfcp_prep.sh
```

`biallelic_snvs.sh` deliberately does not apply a minor-allele-count or
frequency threshold. `split_ref_target.py` applies the minimum minor allele
count of one only after constructing the reference panel. The complete target
VCF in `ref_target/` remains available as truth;
`mask_target/target_observed.vcf.gz` contains the array-observed markers and
`target_masked_true.vcf.gz` contains the withheld markers used for evaluation.
