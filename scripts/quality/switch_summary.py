#!/usr/bin/env python3
"""Summarize retained merged switch intervals against complete haplotype lengths."""
import argparse
import re
from collections import defaultdict

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("assembly_fai")
parser.add_argument("switch_bed")
args = parser.parse_args()
lengths, denominators = {}, defaultdict(int)
with open(args.assembly_fai, encoding="utf-8") as handle:
    for line in handle:
        fields = line.split()
        match = re.fullmatch(r"chr[1-8]\.([1-4])", fields[0])
        if match:
            lengths[fields[0]] = int(fields[1])
            denominators[match.group(1)] += int(fields[1])
errors = defaultdict(int)
intervals = defaultdict(list)
with open(args.switch_bed, encoding="utf-8") as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end, *_ = line.split()
        if chrom not in lengths:
            raise ValueError("Switch interval is outside the nuclear haplotype set: " + chrom)
        start, end = int(start), int(end)
        if start < 0 or end <= start:
            raise ValueError("Invalid switch interval")
        intervals[chrom].append((start, end))
for chrom, spans in intervals.items():
    previous_end = -1
    for start, end in sorted(spans):
        if start < previous_end:
            raise ValueError("Supply merged, nonoverlapping switch intervals")
        previous_end = end
        errors[chrom.rsplit('.', 1)[1]] += end - start
print("haplotype\tswitch_region_length_bp\thaplotype_length_bp\tswitch_error_rate")
for h in sorted(denominators):
    denominator = denominators[h]
    print(f"hap{h}\t{errors[h]}\t{denominator}\t{errors[h]/denominator:.10g}")
total = sum(denominators.values())
if total:
    print(f"all\t{sum(errors.values())}\t{total}\t{sum(errors.values())/total:.10g}")
