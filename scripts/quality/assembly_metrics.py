#!/usr/bin/env python3
"""Basic sequence/contig metrics; gaps are maximal runs of one or more Ns."""
import argparse
import json
import re
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "progeny"))
from localize_markers import fasta


def n50(lengths):
    total = sum(lengths)
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative * 2 >= total:
            return length
    return 0


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("fasta")
parser.add_argument("--telomere-ends", help="TSV: chromosome left|right 0|1; one reviewed row per end")
args = parser.parse_args()
sequence_lengths, contig_lengths = [], []
gaps, gap_bases = 0, 0
for name, sequence in fasta(args.fasta):
    sequence_lengths.append(len(sequence))
    runs = list(re.finditer(r"N+", sequence))
    gaps += len(runs)
    gap_bases += sum(m.end() - m.start() for m in runs)
    contig_lengths.extend(len(s) for s in re.split(r"N+", sequence) if s)
result = dict(sequence_count=len(sequence_lengths), total_length_bp=sum(sequence_lengths),
              sequence_N50_bp=n50(sequence_lengths), contig_count=len(contig_lengths),
              contig_N50_bp=n50(contig_lengths), gap_count=gaps, gap_bases=gap_bases)
if args.telomere_ends:
    ends = {}
    with open(args.telomere_ends, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, side, present = line.split()
            if side not in ("left", "right") or present not in ("0", "1"):
                raise ValueError("Expected chromosome left|right 0|1")
            if (chrom, side) in ends:
                raise ValueError("Duplicate chromosome end")
            ends[chrom, side] = int(present)
    result["reviewed_ends"] = len(ends)
    result["telomere_positive_ends"] = sum(ends.values())
print(json.dumps(result, indent=2))
