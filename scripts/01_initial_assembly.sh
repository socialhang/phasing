#!/usr/bin/env bash
# Initial assembly command transcribed from the project's work.sh.
# Software: hifiasm 0.25.0-r726.
# Run within a PBS job or on an allocated compute node (128 threads).
# Usage: HIFIASM=/absolute/path/to/hifiasm bash 01_initial_assembly.sh /path/to/inputs
# All analysis options and the input order below are preserved from the source.

set -euo pipefail

HIFIASM="${HIFIASM:-hifiasm}"
WORK_DIR="${1:-.}"
cd -- "$WORK_DIR"

exec "$HIFIASM" \
  -t 128 \
  -l 0 \
  --hom-cov 340 \
  --n-hap 4 \
  --h1 SRR9026577_1.fq.f.gz,SRR9026578_1.fq.f.gz \
  --h2 SRR9026577_2.fq.f.gz,SRR9026578_2.fq.f.gz \
  --ont \
  round2.fq.gz \
  ul10.pass.fq.gz \
  ul1.pass.fq.gz \
  ul2.pass.fq.gz \
  ul3.pass.fq.gz \
  ul4.pass.fq.gz \
  ul5.pass.fq.gz \
  ul6.pass.fq.gz \
  ul7.pass.fq.gz \
  ul8.pass.fq.gz \
  ul9.pass.fq.gz \
  all.filtered.fq.gz \
  --dual-scaf \
  --telo-m CCCTAAA
