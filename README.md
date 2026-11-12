# Phasing methods

This repository contains scripts and workflow documentation for genome assembly, phasing using progeny k-mers and Pore-C data, and downstream pangenome and association analyses in autotetraploid alfalfa.

## Workflow manual

The integrated [Supplementary pipeline manual](docs/supplementary-pipeline-manual.md) covers the nine assembly and phasing stages, with an overview diagram, execution parameters, commands, and input/output conventions. An [editable Word version](docs/supplementary-pipeline-manual.docx) is available for manuscript submission.

| Step | Documentation | Script |
| --- | --- | --- |
| 1. Initial assembly | [hifiasm version, inputs, command, and parameters](docs/01-initial-assembly.md) | [01_initial_assembly.sh](scripts/01_initial_assembly.sh) |
| 2. Graph-node depth evaluation | [UL depth estimates, collapsed-node identification, and copy multiplicity](docs/02-node-depth.md) | [02_node_depth.sh](scripts/02_node_depth.sh) |
| 3. Manual chromosome-group partitioning | [Reference alignment, Bandage inspection, and segregation rationale](docs/03-chromosome-groups.md) | [Reference alignment](scripts/03_reference_alignment.sh); manual grouping in Bandage |
| 4. Progeny marker processing and S1 anchoring | [51-mer inheritance vectors and 3D-DNA anchoring](docs/04-progeny-phasing.md) | [04_progeny_kmers.sh](scripts/04_progeny_kmers.sh) |
| 5. Pore-C data analysis | [C-phasing mapping, partitioning, scaffolding, and heatmaps](docs/05-porec-phasing.md) | [05_porec_phasing.sh](scripts/05_porec_phasing.sh) |
| 6. Graph resolution | [Collapsed-node duplication, edge pruning, artifact removal, and sequence merging](docs/06-graph-resolution.md) | Manual operations in Bandage |
| 7. Haplotype-specific gap filling | [Unique 51-mer anchors, UL reads, and four gap closures](docs/07-gap-filling.md) | [07_gap_filling.sh](scripts/07_gap_filling.sh) and reviewed fragment plans |
| 8. 45S rDNA array construction | [Consensus repeat assembly, copy-number evaluation, and uniform Chr6 allocation](docs/08-rdna-assembly.md) | [08_rdna_assembly.sh](scripts/08_rdna_assembly.sh) |
| 9. Genome quality assessment | [Assembly statistics, read support, completeness, QV, and switch errors](docs/09-genome-quality.md) | [09_genome_quality.sh](scripts/09_genome_quality.sh) |
| 10. Pangenome and genome-wide association analyses | [Pangenome construction, variant QC, dosage-resolved GWAS and candidate gene prioritization](docs/10-pangenome-construction.md) | Command examples in the guide |

## System requirements

- Initial assembly: hifiasm **0.25.0-r726**, Bash, and a compute-node allocation suitable for 128 threads. See the [assembly guide](docs/01-initial-assembly.md).
- Node depth: GraphAligner, bedtools, Python 3 with [the helper dependencies](requirements-node-depth.txt); R with ggplot2 and ggpmisc for optional plots. See the [node-depth guide](docs/02-node-depth.md).
- Chromosome grouping: minimap2 for reference alignment and Bandage **v0.9.0** for manual graph inspection; Python 3 for the optional node-list helper.
- Progeny phasing: FASTK **v1.2**, 3D-DNA **201008**, NumPy/SciPy, and the tools listed in the [Step 4 guide](docs/04-progeny-phasing.md).
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
The `03.filter_vcf.py` example reports site quality, depth and genotype counts for the bundled `chr7.1` VCF fragment. Pass the maximum missingness fraction as the second argument. See the [variant calling and QC guide](docs/10-pangenome-construction.md#b-variant-calling-and-quality-control) for the full-chromosome workflow.

```bash
python3 03.filter_vcf.py Demos/raw.vcf.gz 0.1 > demo.filtered.vcf
```
