#!/usr/bin/env bash
set -euo pipefail

data_dir='data/1000g_phase3_v5b'
windows_dir="${1:-$data_dir/windows}"

to_bitpacked() {
    local vcf=$1
    local output=$2
    local n_samples
    local n_haplotypes
    local n_loci
    local tmp
    n_samples=$(bcftools query -l "$vcf" | wc -l)
    n_haplotypes=$((2 * n_samples))
    n_loci=$(bcftools index --nrecords "$vcf")
    tmp=$(mktemp "${output}.XXXXXX")

    echo "bitpacking: $vcf"
    # VCF rows are loci. Write each as ceil(N/64) little-endian words; the
    # SeqArray reader transposes them into sequence-major storage while loading.
    if ! bcftools query -f '[%GT]\n' "$vcf" \
        | perl -e '
            use strict;
            use warnings;

            my ($n, $l) = @ARGV;
            print pack("a4l<l<", "DFCP", $n, $l);

            my $padding = "0" x ((64 - $n % 64) % 64);
            my $loci = 0;
            while (my $line = <STDIN>) {
                $line =~ s/[|\s]//g;
                die "invalid non-binary allele at locus $loci\n" if $line =~ /[^01]/;
                die "locus $loci has " . length($line) . " haplotypes; expected $n\n"
                    if length($line) != $n;

                print pack("b*", $line . $padding);
                ++$loci;
            }

            die "read $loci loci; expected $l\n" if $loci != $l;
        ' "$n_haplotypes" "$n_loci" > "$tmp"
    then
        rm -f -- "$tmp"
        return 1
    fi

    mv -- "$tmp" "$output"
}

write_observed_loci() {
    local ref_vcf=$1
    local target_observed_vcf=$2
    local output=$3
    local n_target_loci
    local target_variants
    local tmp
    n_target_loci=$(bcftools index --nrecords "$target_observed_vcf")
    target_variants=$(mktemp "${output}.target.XXXXXX")
    tmp=$(mktemp "${output}.XXXXXX")

    if ! bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$target_observed_vcf" > "$target_variants"
    then
        rm -f -- "$target_variants" "$tmp"
        return 1
    fi

    if ! bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$ref_vcf" \
        | awk -v expected="$n_target_loci" '
            NR == FNR {
                if ($0 in target_idx) {
                    print "duplicate target variant: " $0 > "/dev/stderr"
                    exit 1
                }
                target_idx[$0] = FNR - 1
                next
            }
            $0 in target_idx {
                if (target_idx[$0] != found) {
                    print "target variants are not in reference order" > "/dev/stderr"
                    exit 1
                }
                print FNR - 1
                delete target_idx[$0]
                ++found
            }
            END {
                if (found != expected) {
                    print "matched " found " target variants; expected " expected > "/dev/stderr"
                    exit 1
                }
            }
        ' "$target_variants" - > "$tmp"
    then
        rm -f -- "$target_variants" "$tmp"
        return 1
    fi

    rm -f -- "$target_variants"
    mv -- "$tmp" "$output"
}

shopt -s nullglob
window_dirs=("$windows_dir"/window_*)
if (( ${#window_dirs[@]} == 0 )); then
    echo "error: no window directories in $windows_dir" >&2
    exit 1
fi

for window_dir in "${window_dirs[@]}"; do
    ref_vcf="$window_dir/ref.vcf.gz"
    target_observed_vcf="$window_dir/target_observed.vcf.gz"
    target_masked_true_vcf="$window_dir/target_masked_true.vcf.gz"
    for vcf in "$ref_vcf" "$target_observed_vcf" "$target_masked_true_vcf"; do
        if [[ ! -f "$vcf" ]]; then
            echo "error: missing window VCF: $vcf" >&2
            exit 1
        fi
    done

    to_bitpacked "$ref_vcf" "$window_dir/ref.bin"
    to_bitpacked "$target_observed_vcf" "$window_dir/target_observed.bin"
    to_bitpacked "$target_masked_true_vcf" "$window_dir/target_masked_true.bin"
    write_observed_loci "$ref_vcf" "$target_observed_vcf" "$window_dir/observed_loci.txt"
done
