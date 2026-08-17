#!/usr/bin/env bash
set -euo pipefail

data_dir='data/1000g_phase3_v5b'
out_dir="$data_dir/dfcp_prep"

mkdir -p "$out_dir"

to_dense_transpose() {
    local vcf=$1
    local name
    name=$(basename "${vcf%.vcf.gz}")

    echo "densifying and transposing: $vcf"
    bcftools query -f '[%GT ]\n' "$vcf" \
        | tr '|' ' ' \
        | sed 's/ $//' \
        | datamash transpose -t ' ' \
        | tr -d ' ' \
        > "$out_dir/$name.txt"
}

to_dense_transpose "$data_dir/mask_target/target_masked_true.vcf.gz"
to_dense_transpose "$data_dir/mask_target/target_masked.vcf.gz"
to_dense_transpose "$data_dir/ref_target/ref.vcf.gz"
