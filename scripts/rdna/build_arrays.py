#!/usr/bin/env python3
"""Build four equal tandem arrays from the 10,261 bp 45S consensus unit.

The assembly-construction default is 781 copies per haplotype (3,124 total).
"""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "progeny"))
from localize_markers import fasta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("unit_fasta")
    parser.add_argument("copies_per_haplotype", type=int, nargs="?", default=781,
                        help="Copies per Chr6 haplotype (default: 781; 3,124 total)")
    parser.add_argument("--fasta", default="rdna_arrays.fa")
    parser.add_argument("--summary", default="rdna_arrays.tsv")
    args = parser.parse_args()
    records = list(fasta(args.unit_fasta))
    if len(records) != 1:
        parser.error("Supply exactly one consensus repeat unit")
    unit_id, sequence = records[0]
    if len(sequence) != 10261 or set(sequence) - set("ACGT"):
        parser.error("The consensus must contain exactly 10,261 A/C/G/T bases")
    if args.copies_per_haplotype < 1:
        parser.error("Copy count must be a positive integer")
    if Path(args.fasta).resolve() == Path(args.summary).resolve():
        parser.error("FASTA and summary paths must differ")
    if Path(args.unit_fasta).resolve() in (Path(args.fasta).resolve(), Path(args.summary).resolve()):
        parser.error("Output paths must differ from the input FASTA")
    with open(args.fasta, "w", encoding="utf-8") as output, open(args.summary, "w", encoding="utf-8") as summary:
        summary.write("chromosome\tarray_id\tunit_id\tunit_length_bp\tcopies\tarray_length_bp\n")
        for haplotype in range(1, 5):
            chromosome = f"chr6.{haplotype}"
            array_id = f"{chromosome}_45S_consensus_array"
            output.write(f">{array_id}\n")
            for _ in range(args.copies_per_haplotype):
                output.write(sequence + "\n")
            length = len(sequence) * args.copies_per_haplotype
            summary.write(f"{chromosome}\t{array_id}\t{unit_id}\t{len(sequence)}\t{args.copies_per_haplotype}\t{length}\n")


if __name__ == "__main__":
    main()
