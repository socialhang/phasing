#!/usr/bin/env python3
"""Select position-labelled unique 51-mers fully contained in flank BED intervals.

Usage: select_flank_markers.py flanks.bed genome_unique_51mers.fa
FASTA identifiers: sequence_id_0basedStart_strand, as in Step 4.
"""
import bisect
import sys

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "progeny"))
from localize_markers import fasta


def main():
    intervals = {}
    with open(sys.argv[1], encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.split()
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError("Invalid BED interval: " + line.strip())
            intervals.setdefault(chrom, []).append((start, end))
    merged = {}
    for chrom, spans in intervals.items():
        union = []
        for start, end in sorted(spans):
            if union and start <= union[-1][1]:
                union[-1] = (union[-1][0], max(end, union[-1][1]))
            else:
                union.append((start, end))
        merged[chrom] = ([start for start, _ in union], union)
    for name, sequence in fasta(sys.argv[2]):
        chrom, start, strand = name.rsplit("_", 2)
        if len(sequence) != 51:
            raise ValueError("Expected a 51-mer: " + name)
        if chrom not in merged:
            continue
        start = int(start)
        starts, spans = merged[chrom]
        index = bisect.bisect_right(starts, start) - 1
        if index >= 0 and start + 51 <= spans[index][1]:
            print(f">{name}\n{sequence}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: select_flank_markers.py flanks.bed genome_unique_51mers.fa")
    main()
