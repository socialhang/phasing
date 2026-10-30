#!/usr/bin/env bash
# Run from cphasing_output/2.prepare/
set -euo pipefail

cphasing prepare ../xjdy_old.fa ../porec_all.pairs.gz -p None -q 0 -t 100 --skip-pairs2contacts
