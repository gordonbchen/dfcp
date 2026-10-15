#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
data_dir="$repo_root/data/1000g_phase3_v5b"
source_dir="$data_dir/ref_target"
mask_dir="$data_dir/mask_target"
filtered_dir="${FILTERED_DIR:-$data_dir/maf1}"
windows_dir="${WINDOWS_DIR:-$data_dir/windows_8_maf1}"
n_windows="${N_WINDOWS:-8}"
overlap="${OVERLAP:-0}"
threads="${THREADS:-1}"

cd "$repo_root"

for input in "$source_dir/ref.vcf.gz" "$source_dir/target.vcf.gz" \
             "$mask_dir/HumanOmni2-5-8-v1-2-A-chr20.bed"; do
    if [[ ! -f "$input" ]]; then
        echo "error: missing input: $input" >&2
        exit 1
    fi
done
if [[ -e "$windows_dir" ]] && find "$windows_dir" -mindepth 1 -print -quit | grep -q .; then
    echo "error: windows output directory must be empty: $windows_dir" >&2
    exit 1
fi
mkdir -p "$filtered_dir"

ref="$filtered_dir/ref.vcf.gz"
target="$filtered_dir/target.vcf.gz"
observed="$filtered_dir/target_observed.vcf.gz"
truth="$filtered_dir/target_masked_true.vcf.gz"
bed="$mask_dir/HumanOmni2-5-8-v1-2-A-chr20.bed"

if [[ -f "$ref" && -f "$target" && -f "$observed" && -f "$truth" ]]; then
    echo "Reusing filtered VCFs in $filtered_dir."
elif find "$filtered_dir" -mindepth 1 -print -quit | grep -q .; then
    echo "error: incomplete filtered VCF set in $filtered_dir" >&2
    exit 1
else
    echo "Filtering reference loci to MAF >= 1% among reference samples."
    bcftools view --threads "$threads" -q 0.01:minor -Oz -o "$ref" "$source_dir/ref.vcf.gz"
    bcftools index --force --threads "$threads" "$ref"

    # Keep target genotypes only at the selected reference loci. Filtering both
    # panels independently would use different allele frequencies and break truth
    # alignment.
    echo "Intersecting target genotypes with filtered reference loci."
    bcftools isec --threads "$threads" -n=2 -w1 -Oz -o "$target" \
        "$source_dir/target.vcf.gz" "$ref"
    bcftools index --force --threads "$threads" "$target"

    echo "Splitting target loci into observed markers and masked truth."
    bcftools view --threads "$threads" -T "$bed" -Oz -o "$observed" "$target"
    bcftools index --force --threads "$threads" "$observed"
    bcftools view --threads "$threads" -T ^"$bed" -Oz -o "$truth" "$target"
    bcftools index --force --threads "$threads" "$truth"
fi

ref_count=$(bcftools index --nrecords "$ref")
target_count=$(bcftools index --nrecords "$target")
observed_count=$(bcftools index --nrecords "$observed")
truth_count=$(bcftools index --nrecords "$truth")
expected_ref_count=$(bcftools query -f '%AC\t%AN\n' "$source_dir/ref.vcf.gz" | awk -F '\t' '
    {
        mac = ($1 < $2 - $1) ? $1 : $2 - $1
        if (100 * mac >= $2) ++n
    }
    END { print n + 0 }
')
if (( ref_count != expected_ref_count )); then
    echo "error: MAF filter retained $ref_count loci; expected $expected_ref_count from reference genotypes" >&2
    exit 1
fi
if (( ref_count == 0 || ref_count != target_count || ref_count != observed_count + truth_count )); then
    echo "error: filtered locus counts do not align: ref=$ref_count target=$target_count observed=$observed_count truth=$truth_count" >&2
    exit 1
fi

# Check the actual called reference genotypes, with integer arithmetic so the
# 1% boundary is exact for the available allele counts.
bcftools query -f '%AC\t%AN\n' "$ref" | awk -F '\t' '
    NF != 2 || $1 !~ /^[0-9]+$/ || $2 !~ /^[0-9]+$/ { exit 2 }
    {
        mac = ($1 < $2 - $1) ? $1 : $2 - $1
        if (100 * mac < $2) {
            print "error: retained site below 1% MAF (MAC=" mac ", AN=" $2 ")" > "/dev/stderr"
            exit 1
        }
    }
' || { echo "error: failed MAF validation" >&2; exit 1; }

mkdir -p "$windows_dir"
python3 "$script_dir/window.py" \
    --ref "$ref" \
    --target-observed "$observed" \
    --target-masked-true "$truth" \
    --output-dir "$windows_dir" \
    --n-windows "$n_windows" \
    --overlap "$overlap" \
    --threads "$threads"
"$script_dir/dfcp_prep.sh" "$windows_dir"

printf 'Prepared MAF>=1%% panel: %s loci (observed=%s, masked=%s), %s windows at %s\n' \
    "$ref_count" "$observed_count" "$truth_count" "$n_windows" "$windows_dir"
