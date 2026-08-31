# 1000 Genomes Phase 3 v5b preprocessing

The pipeline keeps VCF as the inspectable source format through windowing. DFCP binaries are generated
only after each reference and target window has been written and checked.

## Complete pipeline

Run from the repository root:

```bash
THREADS=8 N_WINDOWS=200 OVERLAP=32 N_GENERATE=1 \
  WINDOWS_DIR=data/1000g_phase3_v5b/windows_200_o32 \
  scripts/1000g_phase3_v5b/prep_data.sh
```

`N_WINDOWS` and `OVERLAP` define the chromosome-wide window plan. `N_GENERATE` writes only the first
requested windows, which is useful during development; it defaults to all `N_WINDOWS`. The output
directory must not exist or must be empty.

The stages are:

1. `download.sh` downloads the chr20 VCF, sample panel, Omni array BED, and genetic map.
2. `biallelic_snvs.sh` retains PASS records with one-base REF and ALT alleles.
3. `split_ref_target.py` selects target individuals, recalculates reference-panel `AC` and `AN`, removes
   reference-monomorphic records, and aligns the complete target VCF to the reference records.
4. `mask_target.sh` partitions target records into observed array markers and masked truth.
5. `window.py` writes aligned overlapping VCF windows.
6. `dfcp_prep.sh` bitpacks every generated window and writes its local observed-locus indexes.

## Window generation

Run the VCF window stage directly with:

```bash
python3 scripts/1000g_phase3_v5b/window.py \
  --output-dir data/1000g_phase3_v5b/windows_200_o32 \
  --n-windows 200 --overlap 32 --n-generate 1 --threads 8
```

Every generated `window_NNNN/` initially contains:

```text
ref.vcf.gz
target_observed.vcf.gz
target_masked_true.vcf.gz
```

Each VCF is indexed. `windows.tsv` describes every planned window and marks which ones were generated.
Window sizes are planned from reference-record indexes, then each window is extracted as the physical
interval between its first and last record. The same interval is applied to all three VCFs.

Bitpack the generated windows with:

```bash
scripts/1000g_phase3_v5b/dfcp_prep.sh data/1000g_phase3_v5b/windows_200_o32
```

This adds `ref.bin`, `target_observed.bin`, `target_masked_true.bin`, and `observed_loci.txt` to each
window directory. The observed-loci file contains zero-based indexes into that window's reference VCF.

## Window visualization

Reproduce the interactive selection report directly from the chromosome-wide reference and observed
target VCFs:

```bash
.venv/bin/python scripts/1000g_phase3_v5b/window_viz.py --output docs/window.html
```

Use `--physical-only` to build the report in Mb without the genetic map.

## Imputation and evaluation

```bash
scripts/1000g_phase3_v5b/impute.sh data/1000g_phase3_v5b/windows_200_o32 \
  --n-parallel 4 --mode soft --init pbwt --viterbi_impute 0
./build/eval_impute data/1000g_phase3_v5b/windows_200_o32 output/imputation.tsv
.venv/bin/python scripts/impute_viz.py output/imputation.tsv --output output/imputation.html
```

`impute.sh` builds DFCP once, runs the same options on every window, and writes `probs.bin`, `impute.json`,
and `impute.log` in each window directory. `--n-parallel` controls the number of simultaneous DFCP processes
and defaults to one. Each process's stderr is printed as one block when it finishes. `Ctrl-C` stops every
active process. `eval_impute` reads `AC` and `AN` from each reference VCF, keeps each overlap locus from the
window where it is farthest from an edge, and pools r-squared and accuracy over all retained target alleles
at each MAC.
