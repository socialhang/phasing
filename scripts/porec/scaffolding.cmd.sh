#!/usr/bin/env bash
# Run from cphasing_output/4.scaffolding/
set -euo pipefail

cphasing scaffolding \
    ../3.hyperpartition/output.clusters.txt ../2.prepare/porec_all.counts_AAGCTT.txt \
     ../2.prepare/porec_all.clm.gz   \
     -at ../1.alleles/xjdy_old.allele.table \
     -sc ../2.prepare/porec_all.split.contacts \
     -f ../xjdy_old.fa \
     -t 100 \
     -o groups.agp \
     -m precision
