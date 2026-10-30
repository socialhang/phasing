#!/usr/bin/env python3
"""Localize canonical unique 51-mers; emit node_0basedStart_strand FASTA.

Portable sequence-to-coordinate helper supplied with the workflow manual.
"""
import sys


def fasta(path):
    name, parts = None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: localize_markers.py parents.fa unique_51mers.fa")
    complement = str.maketrans("ACGT", "TGCA")
    markers = {seq for _, seq in fasta(sys.argv[2])}
    if any(len(seq) != 51 or set(seq) - set("ACGT") for seq in markers):
        raise SystemExit("Marker sequences must be canonical A/C/G/T 51-mers")
    seen_ids = set()
    for name, sequence in fasta(sys.argv[1]):
        if name in seen_ids:
            raise SystemExit("Duplicate parental sequence ID: " + name)
        seen_ids.add(name)
        for start in range(len(sequence) - 50):
            forward = sequence[start:start + 51]
            if set(forward) - set("ACGT"):
                continue
            reverse = forward.translate(complement)[::-1]
            canonical = min(forward, reverse)
            if canonical in markers:
                strand = "+" if canonical == forward else "-"
                print(f">{name}_{start}_{strand}\n{canonical}")


if __name__ == "__main__":
    main()
