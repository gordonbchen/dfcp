#!/usr/bin/env bash
set -euo pipefail

data_dir='data/1000g_phase3_v5b'
output_dir="$data_dir/biallelic_snvs"
default_input="$data_dir/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
default_output="$output_dir/biallelic_snvs.vcf.gz"

input_vcf="${1:-$default_input}"
output_vcf="${2:-$default_output}"
threads="${THREADS:-1}"
sv_report="${output_vcf%.vcf.gz}.svtype_counts.tsv"

if ! command -v bcftools >/dev/null; then
    echo 'error: bcftools is required' >&2
    exit 1
fi
if [[ ! -f "$input_vcf" ]]; then
    echo "error: input VCF does not exist: $input_vcf" >&2
    exit 1
fi
if [[ "$output_vcf" != *.vcf.gz ]]; then
    echo 'error: output must end in .vcf.gz' >&2
    exit 1
fi

# count structural variant types.
mkdir -p "$(dirname "$output_vcf")"

{
    printf 'svtype\tcount\n'
    bcftools query -f '%INFO/SVTYPE\n' "$input_vcf" |
        awk '$1 != "." { counts[$1]++ } END { for (type in counts) print type "\t" counts[type] }' |
        sort
} > "$sv_report"

# -f PASS: retain records that passed the call-set filters.
# -m2 -M2: require exactly two alleles (one REF and one ALT).
# -v snps: require both alleles to be one base long; this excludes indels and SVs.
bcftools view \
    --threads "$threads" \
    -f PASS \
    -m2 -M2 \
    -v snps \
    -Oz \
    -o "$output_vcf" \
    "$input_vcf"
bcftools index --force "$output_vcf"

input_records=$(bcftools query -f '%CHROM\n' "$input_vcf" | wc -l)
output_records=$(bcftools query -f '%CHROM\n' "$output_vcf" | wc -l)

printf 'Input records: %s\n' "$input_records"
printf 'Retained PASS biallelic SNVs: %s\n' "$output_records"
printf 'Filtered VCF: %s\n' "$output_vcf"
printf 'Structural-variant audit: %s\n' "$sv_report"
