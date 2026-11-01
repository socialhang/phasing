#!/usr/bin/env python3
"""Normalize Jellyfish marker counts against a dataset-specific baseline.

Input: two whitespace-separated columns, k-mer and count (jellyfish query output).
"""
import argparse
import math


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts")
    parser.add_argument("baseline", type=float)
    parser.add_argument("--baseline-copies", type=float, default=1.0,
                        help="Copies represented by the normalization markers in the tetraploid genome")
    args = parser.parse_args()
    if not math.isfinite(args.baseline) or args.baseline <= 0:
        parser.error("Baseline must be finite and positive")
    if not math.isfinite(args.baseline_copies) or args.baseline_copies <= 0:
        parser.error("Baseline copy multiplicity must be finite and positive")
    print("kmer\tcount\tdepth_ratio\tscaled_copy_abundance")
    with open(args.counts, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            marker, count = line.split()
            count = int(count)
            if count < 0:
                raise ValueError("Negative marker count")
            ratio = count / args.baseline
            print(f"{marker}\t{count}\t{ratio:.8g}\t{ratio * args.baseline_copies:.8g}")


if __name__ == "__main__":
    main()
