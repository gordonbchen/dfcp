#!/usr/bin/env python3

import argparse
import os
import struct
import tempfile
from pathlib import Path


SEQ_HEADER = struct.Struct("<4sII")
WORD = struct.Struct("<Q")


def temporary_path(output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{output.name}.", dir=output.parent)
    os.close(fd)
    return Path(name)


def pack_locus(alleles: list[str], output) -> None:
    for start in range(0, len(alleles), 64):
        word = 0
        for bit, allele in enumerate(alleles[start:start + 64]):
            if allele not in ("0", "1"):
                raise ValueError(f"invalid non-binary allele: {allele}")
            word |= int(allele) << bit
        output.write(WORD.pack(word))


def prepare(gen_file: Path, output_dir: Path) -> tuple[int, int]:
    seq_output = output_dir / "ref.bin"
    pos_output = output_dir / "variant_pos.txt"
    seq_tmp = temporary_path(seq_output)
    pos_tmp = temporary_path(pos_output)

    try:
        with gen_file.open() as source, seq_tmp.open("wb") as sequences:
            header = source.readline().split()
            if header[:4] != ["Chrom", "Pos", "Anc_all", "Der_all"] or len(header) <= 4:
                raise ValueError("invalid fastsimcoal .gen header")

            n_haplotypes = len(header) - 4
            positions = []
            chromosome = None
            sequences.write(SEQ_HEADER.pack(b"DFCP", n_haplotypes, 0))

            for line_number, line in enumerate(source, 2):
                if not line.strip():
                    continue
                fields = line.split()
                if len(fields) != len(header):
                    raise ValueError(f"line {line_number} has {len(fields)} columns; expected {len(header)}")

                if chromosome is None:
                    chromosome = fields[0]
                elif fields[0] != chromosome:
                    raise ValueError("DFCP position files support one chromosome")

                position = int(fields[1])
                if positions and position < positions[-1]:
                    raise ValueError("variant positions are not sorted")
                positions.append(position)
                pack_locus(fields[4:], sequences)

            if not positions:
                raise ValueError("fastsimcoal .gen file has no variants")
            sequences.seek(0)
            sequences.write(SEQ_HEADER.pack(b"DFCP", n_haplotypes, len(positions)))

        with pos_tmp.open("w") as positions_file:
            positions_file.write(f"{len(positions)}\n{', '.join(map(str, positions))}\n")

        os.replace(seq_tmp, seq_output)
        os.replace(pos_tmp, pos_output)
        return n_haplotypes, len(positions)
    finally:
        seq_tmp.unlink(missing_ok=True)
        pos_tmp.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a haploid fastsimcoal .gen file to DFCP sequence and position files."
    )
    parser.add_argument("gen_file", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    n_haplotypes, n_loci = prepare(args.gen_file, args.output_dir)
    print(f"wrote {args.output_dir / 'ref.bin'}: N={n_haplotypes}, L={n_loci}")
    print(f"wrote {args.output_dir / 'variant_pos.txt'}")


if __name__ == "__main__":
    main()
