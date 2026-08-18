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
   panel, and the Illumina Omni2.5 BED into `data/1000g_phase3_v5b/`.
2. `biallelic_snps.sh` retains PASS biallelic SNPs and writes the filtered VCF,
   CSI index, and structural-variant audit into `biallelic_snps/`.
3. `split_ref_target.py` selects two target individuals per population and
   writes the reference and target VCFs, indexes, and sample lists into
   `ref_target/`.
4. `mask_target.sh` writes the chromosome 20 Omni BED, observed target VCF,
   withheld target truth VCF, and indexes into `mask_target/`.
5. `dfcp_prep.sh` streams each VCF directly into locus-major, 64-bit-padded
   DFCP inputs in `dfcp_prep/`. It also writes
   `target_masked.unmasked_loci.txt`, which maps observed target columns to
   reference loci. The C++ reader performs the packed transpose.

Run stages individually with their defaults if needed:

```bash
scripts/1000g_phase3_v5b/download.sh
scripts/1000g_phase3_v5b/biallelic_snps.sh
python3 scripts/1000g_phase3_v5b/split_ref_target.py --seed 0 --threads 8
THREADS=8 scripts/1000g_phase3_v5b/mask_target.sh
scripts/1000g_phase3_v5b/dfcp_prep.sh
```

`biallelic_snps.sh` deliberately does not apply a minor-allele-count or
frequency threshold. `split_ref_target.py` applies the minimum minor allele
count of one only after constructing the reference panel. The complete target
VCF in `ref_target/` remains available as truth; `mask_target/target_masked.vcf.gz`
contains the array-observed markers and `target_masked_true.vcf.gz` contains the
withheld markers used for evaluation.
