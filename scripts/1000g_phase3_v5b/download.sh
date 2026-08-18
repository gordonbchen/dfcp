#!/usr/bin/env bash
set -euo pipefail

data_dir='data/1000g_phase3_v5b'
force=0
if [[ ${1:-} == '--force' ]]; then
    force=1
    shift
fi
if (( $# != 0 )); then
    echo "usage: $0 [--force]" >&2
    exit 2
fi

mkdir -p "$data_dir"

download() {
    local url=$1
    local destination=$2
    local partial="$destination.part"

    if [[ -e "$destination" && $force -eq 0 ]]; then
        echo "Skipping existing file: $destination"
        return
    fi

    echo "Downloading $url"
    curl --fail --location --retry 3 --output "$partial" "$url"
    mv -- "$partial" "$destination"
}

download \
    'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz' \
    "$data_dir/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

download \
    'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi' \
    "$data_dir/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"

download \
    'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel' \
    "$data_dir/integrated_call_samples_v3.20130502.ALL.panel"

download \
    'https://webdata.illumina.com/downloads/productfiles/humanomni25/v1-2/SupportFiles/HumanOmni2-5-8-v1-2-A.bed' \
    "$data_dir/HumanOmni2-5-8-v1-2-A.bed"

# Download GRCh37 map.
map_archive="$data_dir/plink.GRCh37.map.zip"
map_file="$data_dir/plink.chr20.GRCh37.map"
download \
    'https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip' \
    "$map_archive"

expected_map_sha256='8e88141345844571bc63efac6c043ee29e8ca41f26087dc51b7f1918943a077b'
actual_map_sha256=$(sha256sum "$map_archive" | cut -d ' ' -f 1)
if [[ $actual_map_sha256 != "$expected_map_sha256" ]]; then
    echo "Unexpected SHA-256 for $map_archive" >&2
    exit 1
fi

if [[ ! -e "$map_file" || $force -eq 1 ]]; then
    map_partial="$map_file.part"
    echo "Extracting chr20 GRCh37 genetic map"
    unzip -p "$map_archive" plink.chr20.GRCh37.map > "$map_partial"
    mv -- "$map_partial" "$map_file"
else
    echo "Skipping existing file: $map_file"
fi

echo "Downloads available in $data_dir"
