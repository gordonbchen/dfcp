#!/usr/bin/env bash
set -euo pipefail

data_dir='data/1000g_phase3_v5b'
out_dir="$data_dir/dfcp_prep"

mkdir -p "$out_dir"

to_bitpacked() {
    local vcf=$1
    local name
    local n_samples
    local n_haplotypes
    local n_loci
    local output
    local tmp
    name=$(basename "${vcf%.vcf.gz}")
    n_samples=$(bcftools query -l "$vcf" | wc -l)
    n_haplotypes=$((2 * n_samples))
    n_loci=$(bcftools index --nrecords "$vcf")
    output="$out_dir/$name.bin"
    tmp=$(mktemp "$out_dir/.${name}.XXXXXX")

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
    local n_target_loci
    local target_positions
    local output
    local tmp
    n_target_loci=$(bcftools index --nrecords "$target_observed_vcf")
    target_positions=$(mktemp "$out_dir/.target_positions.XXXXXX")
    output="$out_dir/observed_loci.txt"
    tmp=$(mktemp "$out_dir/.observed_loci.XXXXXX")

    if ! bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$target_observed_vcf" > "$target_positions"
    then
        rm -f -- "$target_positions" "$tmp"
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
        ' "$target_positions" - > "$tmp"
    then
        rm -f -- "$target_positions" "$tmp"
        return 1
    fi

    rm -f -- "$target_positions"
    mv -- "$tmp" "$output"
}

ref_vcf="$data_dir/ref_target/ref.vcf.gz"
target_observed_vcf="$data_dir/mask_target/target_observed.vcf.gz"
target_masked_true_vcf="$data_dir/mask_target/target_masked_true.vcf.gz"

to_bitpacked "$ref_vcf"
to_bitpacked "$target_observed_vcf"
to_bitpacked "$target_masked_true_vcf"
write_observed_loci "$ref_vcf" "$target_observed_vcf"
