#!/usr/bin/env bash
# Run from cphasing_output/5.plot/
set -euo pipefail

min_quality=1
cool_binsize=10k
heatmap_binsize=500k


cphasing pairs2cool ../porec_all.pairs.gz ../xjdy_old.contigsizes porec_all.q${min_quality}.${cool_binsize}.cool -q ${min_quality} -bs ${cool_binsize}
cphasing plot -a ../4.scaffolding/groups.agp -m porec_all.q${min_quality}.${cool_binsize}.cool -o groups.q${min_quality}.${heatmap_binsize}.wg.png -bs ${heatmap_binsize} -oc
