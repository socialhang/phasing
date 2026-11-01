#!/usr/bin/env python3
"""Write selected markers and their reverse complements, one sequence per line."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "progeny"))
from localize_markers import fasta

complement = str.maketrans("ACGT", "TGCA")
for _, seq in fasta(sys.argv[1]):
    print(seq)
    reverse = seq.translate(complement)[::-1]
    if reverse != seq:
        print(reverse)
