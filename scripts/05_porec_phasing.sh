#!/usr/bin/env bash
# Submit through PBS with resources suitable for 100 threads.
# Usage: bash 05_porec_phasing.sh ASSEMBLY_FASTA POREC_FASTQ [WORK_DIR]
set -euo pipefail

ASSEMBLY="$(realpath -- "${1:?Specify the assembly FASTA}")"
POREC="$(realpath -- "${2:?Specify the Pore-C FASTQ}")"
WORK_DIR="${3:-.}"
CPHASING="${CPHASING:-cphasing}"
mkdir -p -- "$WORK_DIR"
cd -- "$WORK_DIR"

"$CPHASING" pipeline -f "$ASSEMBLY" -pcd "$POREC" -t 100 -n 8:4
