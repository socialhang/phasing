#!/usr/bin/env bash
# Supporting steps for manually reviewed gap closure. Submit compute jobs via PBS.
# Usage: bash 07_gap_filling.sh {markers|reads|align|assemble} ...
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
MINIMAP2="${MINIMAP2:-minimap2}"
STAGE="${1:?Specify markers, reads, align, or assemble}"
shift
case "$STAGE" in
  markers)
    BED="${1:?Specify the 50 kb flank BED}"
    CATALOG="${2:?Specify coordinate-labelled genome-wide unique 51-mers}"
    "$PYTHON" "$SCRIPT_DIR/gap_filling/select_flank_markers.py" "$BED" "$CATALOG" > flank_markers.fa
    "$PYTHON" "$SCRIPT_DIR/gap_filling/marker_sequences.py" flank_markers.fa > flank_markers.txt
    ;;
  reads)
    IDS="${1:?Specify candidate read IDs}"
    FASTQ="${2:?Specify UL reads in gzip FASTQ}"
    MARKERS="${3:?Specify flank_markers.txt}"
    "$PYTHON" "$SCRIPT_DIR/gap_filling/filter_have_kmer_reads.py" "$IDS" "$FASTQ" "$MARKERS" > gap_reads.fq
    ;;
  align)
    TARGET="${1:?Specify the flank or assembly FASTA}"
    READS="${2:?Specify the selected reads FASTA/FASTQ}"
    "$MINIMAP2" -x map-ont "$TARGET" "$READS" > gap_reads.paf
    ;;
  assemble)
    SOURCES="${1:?Specify FASTA containing unitigs and selected reads}"
    PLAN="${2:?Specify the reviewed fragment BED plan}"
    "$PYTHON" "$SCRIPT_DIR/gap_filling/resolve_bed_coordinates.py" "$PLAN" > fragment.real.bed
    bedtools getfasta -fi "$SOURCES" -bed fragment.real.bed -name -s > selected_fragments.fa
    "$PYTHON" "$SCRIPT_DIR/gap_filling/link_fragment.py" selected_fragments.fa > gap_closed_sequences.fa
    ;;
  *) printf 'Unknown stage: %s\n' "$STAGE" >&2; exit 1 ;;
esac
