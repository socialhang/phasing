#!/usr/bin/env bash
# Supporting alignment for manual chromosome grouping, not automatic partitioning.
# Usage: bash 03_reference_alignment.sh REFERENCE_FASTA UNITIG_FASTA [OUTPUT_PAF]
# Submit through PBS or run on an allocated compute node (40 threads).
set -euo pipefail

REFERENCE="${1:?Specify the reference FASTA}"
UNITIGS="${2:?Specify the unitig FASTA}"
OUTPUT="${3:-scm.paf}"
MINIMAP2="${MINIMAP2:-minimap2}"

"$MINIMAP2" -x asm5 -t 40 "$REFERENCE" "$UNITIGS" > "$OUTPUT"
