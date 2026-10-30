#!/usr/bin/env bash
# Run from cphasing_output/1.alleles/
set -euo pipefail

cphasing alleles -f ../xjdy_old.fa -k 19 -w 19 -m 0.5 -d 0.1 -t 100
