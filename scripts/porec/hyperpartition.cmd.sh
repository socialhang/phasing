#!/usr/bin/env bash
# Run from cphasing_output/3.hyperpartition/
set -euo pipefail

cphasing hyperpartition \
    ../porec_all.porec.gz ../xjdy_old.contigsizes output.clusters.txt --porec \
     --mode phasing -at ../1.alleles/xjdy_old.allele.table \
     -r1 -1.0 -ir1 0.8 \
     -r2 1.0 -ir2 0.8 \
     -as 0.85 -mao 0.3 \
     -mw 0.1 -q1 1 \
     -q2 2 -mc 5 \
     -ml 10000 -ms 500000.0 \
     -t 100 -n 8:4
