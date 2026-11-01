#!/usr/bin/env bash
# Submit assembly/counting stages through PBS with suitable resources.
# Usage: bash 08_rdna_assembly.sh {assemble|count|query|arrays} ...
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
MBG="${MBG:-MBG}"
THREADS="${THREADS:-40}"
STAGE="${1:?Specify assemble, count, query, or arrays}"
shift
case "$STAGE" in
  assemble)
    READS="${1:?Specify rDNA-containing reads}"
    "$MBG" -i "$READS" -o rdna_consensus.gfa -k 4001 -a 4 -u 4 -t "$THREADS"
    ;;
  count)
    OUTPUT="${1:?Specify the output Jellyfish database}"
    shift
    [[ $# -gt 0 ]] || { printf 'Supply at least one uncompressed FASTA/FASTQ input.\n' >&2; exit 1; }
    jellyfish count -m 51 -s 100M -t 40 -o "$OUTPUT" "$@"
    ;;
  query)
    DATABASE="${1:?Specify the Jellyfish database}"
    MARKERS="${2:?Specify a plain-text list of 51-mers}"
    while IFS= read -r marker || [[ -n "$marker" ]]; do
      marker="${marker%$'\r'}"
      [[ -z "$marker" || "$marker" == \#* ]] && continue
      jellyfish query "$DATABASE" "$marker"
    done < "$MARKERS"
    ;;
  arrays)
    UNIT="${1:?Specify the 10261 bp consensus FASTA}"
    COPIES="${2:-781}"
    "$PYTHON" "$SCRIPT_DIR/rdna/build_arrays.py" "$UNIT" "$COPIES" \
      --fasta rdna_arrays.fa --summary rdna_arrays.tsv
    ;;
  *) printf 'Unknown stage: %s\n' "$STAGE" >&2; exit 1 ;;
esac
