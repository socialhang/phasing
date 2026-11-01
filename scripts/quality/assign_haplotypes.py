#!/usr/bin/env python3
"""Retain the first mapper-ranked PAF alignment per read, as in best0 file 1."""
import argparse
import re
from contextlib import ExitStack

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("paf")
parser.add_argument("--prefix", default="hap")
args = parser.parse_args()
seen = set()
with ExitStack() as stack:
    outputs = {h: stack.enter_context(open(f"{args.prefix}{h}.reads.txt", "w", encoding="utf-8"))
               for h in range(1, 5)}
    with open(args.paf, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 12:
                continue
            read, target = fields[0], fields[5]
            if read in seen:
                continue
            seen.add(read)
            match = re.fullmatch(r"chr[1-8]\.([1-4])", target)
            if match:
                outputs[int(match.group(1))].write(read + "\n")
