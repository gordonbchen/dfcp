#!/usr/bin/env bash
set -euo pipefail

data_dir='data/1000g_phase3_v5b'
input_dir="$data_dir/ref_target"
output_dir="$data_dir/mask_target"
threads="${THREADS:-1}"

mkdir -p "$output_dir"

# change chr20 to 20, get start and end pos.
echo 'Extracting chr20 Omni variant pos.'
awk '$1=="chr20" {print 20, $2, $3}' OFS='\t' \
    "$data_dir/HumanOmni2-5-8-v1-2-A.bed" > "$output_dir/HumanOmni2-5-8-v1-2-A-chr20.bed"

# Create target_observed.vcf.gz from array positions.
# -T keep target pos in file.
echo 'Creating observed target.'
bcftools view \
    "$input_dir/target.vcf.gz" \
    -T "$output_dir/HumanOmni2-5-8-v1-2-A-chr20.bed" \
    -Oz -o "$output_dir/target_observed.vcf.gz" \
    --threads "$threads"

bcftools index --force "$output_dir/target_observed.vcf.gz"

# Create target_masked_true.vcf.gz from positions absent from the array.
# -T ^file remove target pos in file.
echo 'Creating masked target truth.'
bcftools view \
    "$input_dir/target.vcf.gz" \
    -T ^"$output_dir/HumanOmni2-5-8-v1-2-A-chr20.bed" \
    -Oz -o "$output_dir/target_masked_true.vcf.gz" \
    --threads "$threads"

bcftools index --force "$output_dir/target_masked_true.vcf.gz"
