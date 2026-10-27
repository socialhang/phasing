# Phasing methods

This repository contains scripts and workflow documentation for genome assembly, progeny k-mer-based phasing, and variant processing in autotetraploid alfalfa.

## Workflow manual

| Step | Documentation | Script |
| --- | --- | --- |
| 1. Initial assembly | [hifiasm version, inputs, command, and parameters](docs/01-initial-assembly.md) | [01_initial_assembly.sh](scripts/01_initial_assembly.sh) |
| 2. Graph-node depth evaluation | [HiFi alignment, coverage, and per-node KDE estimates](docs/02-node-depth.md) | [02_node_depth.sh](scripts/02_node_depth.sh) |
| 3. Manual chromosome-group partitioning | [Bandage inspection and segregation rationale](docs/03-chromosome-groups.md) | Manual operation in Bandage; no script |

## System requirements

- Initial assembly: hifiasm **0.25.0-r726**, Bash, and a compute-node allocation suitable for 128 threads. See the [assembly guide](docs/01-initial-assembly.md).
- Node depth: GraphAligner, bedtools, Python 3 with [the helper dependencies](requirements-node-depth.txt); R with ggplot2 and ggpmisc for optional plots. See the [node-depth guide](docs/02-node-depth.md).
- The Python examples below: Python > 3.7.

## Installation

Download or clone this repository. Install hifiasm separately for the initial assembly step; the repository scripts themselves do not require installation.

## Functions and Usage of the Script
The script 01.merge_kmer_from_NGS.py is used to merge 51-mers from a large number of offspring individuals and generate a presence/absence matrix of these 51-mers across these individuals.
```
python 01.merge_kmer_from_NGS.py Demos/file.list tcg > mergeout.part_tcg.txt
python 01.merge_kmer_from_NGS.py Demos/file.list atg > mergeout.part_atg.txt
python 01.merge_kmer_from_NGS.py Demos/file.list atc > mergeout.part_atc.txt
......
```
The 02.make_pseudo-contact_matrices.py script calculates pairwise Pearson correlation coefficients across positions in the merged presence/absence matrix, and simulates them as pseudo-contact matrices. These matrices are then used in the next step for phasing via the iterative algorithm within the 3D-DNA pipeline.
```
python 02.make_pseudo-contact_matrices.py Demos/data2.txt Demos/data2.txt > out.txt
```
The 03.filter_vcf.py script is used to process the variant calling VCF files of autotetraploids and calculate variant frequencies. The output file is used for subsequent variant filtering.
```
python 03.filter_vcf.py Demos/raw.vcf.gz
```
