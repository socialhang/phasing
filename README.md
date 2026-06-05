# Phasing methods
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
