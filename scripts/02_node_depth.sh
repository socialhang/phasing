#!/usr/bin/env bash
# Run each stage in a PBS job or on an allocated compute node.
# Use separate working directories for HiFi and UL; output names use the hifi example prefix.
# Usage: bash 02_node_depth.sh {map|convert|coverage|kde} WORK_DIR READS_LIST
# The list contains one HiFi or UL FASTQ path per line, relative to WORK_DIR or absolute.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
STAGE="${1:?Specify map, convert, coverage, or kde}"
WORK_DIR="${2:?Specify the working directory}"
MANIFEST="${3:?Specify a HiFi or UL read list}"
GFA="${GFA:-hifiasm.asm.bp.p_utg.gfa}"
FAI="${FAI:-hifiasm.asm.bp.p_utg.fa.fai}"
THREADS="${THREADS:-40}"
PYTHON="${PYTHON:-python3}"
cd -- "$WORK_DIR"

READS=()
while IFS= read -r read_file || [[ -n "$read_file" ]]; do
  read_file="${read_file%$'\r'}"
  [[ -z "$read_file" || "$read_file" == \#* ]] && continue
  READS+=("$read_file")
done < "$MANIFEST"
[[ ${#READS[@]} -gt 0 ]] || { printf 'The read list is empty.\n' >&2; exit 1; }

case "$STAGE" in
  map)
    for read_file in "${READS[@]}"; do
      GraphAligner -t "$THREADS" -g "$GFA" -f "$read_file" \
        -a "$(basename -- "$read_file").gaf" -x vg
    done
    ;;
  convert)
    for read_file in "${READS[@]}"; do
      gaf="$(basename -- "$read_file").gaf"
      "$PYTHON" "$SCRIPT_DIR/node_depth/gaf2paf.py" "$GFA" "$gaf" > "$gaf.bed"
    done
    ;;
  coverage)
    # Explicit inputs keep prior merged BED outputs out of a repeated run.
    for read_file in "${READS[@]}"; do
      cat -- "$(basename -- "$read_file").gaf.bed"
    done > hifi.raw.bed
    awk '{print $3"\t"$4"\t"$5"\t"$1"\t"$6}' hifi.raw.bed \
      | sort -k1,1 -k2,2n > hifi.bed
    bedtools genomecov -i hifi.bed -g "$FAI" -bga > hifi.base.depth
    gzip -c hifi.base.depth > hifi.base.depth.gz
    ;;
  kde)
    "$PYTHON" "$SCRIPT_DIR/node_depth/get_KDE_Mode_v2.py" hifi.base.depth.gz > hifi.last.depth
    ;;
  *)
    printf 'Unknown stage: %s\nChoose map, convert, coverage, or kde.\n' "$STAGE" >&2
    exit 1
    ;;
esac
