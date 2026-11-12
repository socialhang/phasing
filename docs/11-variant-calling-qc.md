# Step 11 Tetraploid variant calling and quality control

Short-read whole-genome sequencing data from 358 autotetraploid alfalfa accessions were aligned to an H1 linear reference. The calling files use one coordinate system for all eight chromosomes (`chr1.1` through `chr8.1`). Each accession was genotyped at ploidy four so the VCF retains tetraploid allele dosage. The GVCF manifest and cohort VCF contain 358 samples. The site-selection script uses 385 as a fixed denominator.

## 1. Align and prepare one BAM per accession

Use a sample manifest linking each accession ID to its paired FASTQ files. Use the same reference FASTA, FASTA index and sequence dictionary for all accessions. The archived calling scripts use `chr.dantao.fa`. Verify the sequence dictionary before joining these positions to another XJDYt2t H1 assembly release.

BWA-MEM aligns the short reads. SAMtools retains alignments with MAPQ at least 1 and excludes secondary (`0x100`) and supplementary (`0x800`) records. The example uses a pair of FASTQ files for one accession and writes a coordinate-sorted BAM.

```bash
bwa mem -t "$threads" "$ref" "$r1" "$r2" \
  | samtools view -u -q 1 -F 0x900 - \
  | samtools sort -@ "$threads" -o "$sample.q1.bam" -
```

Assign the accession ID as the GATK sample name, mark duplicates and index the resulting BAM. The read-group and duplicate-marking options below are from the archived per-accession script.

```bash
gatk AddOrReplaceReadGroups \
  -I "$sample.q1.bam" -O "$sample.rg.bam" \
  -LB library1 -PL illumina -PU pl1 -SM "$sample"
mkdir -p tmp
gatk MarkDuplicates \
  -I "$sample.rg.bam" -O "$sample.markdup.bam" \
  -M "$sample.markdup.metrics.txt" --READ_NAME_REGEX null --TMP_DIR tmp
samtools index -@ "$threads" "$sample.markdup.bam"
```

## 2. Call GVCFs and jointly genotype the cohort

Run HaplotypeCaller in GVCF mode for each accession. The archive uses `-ploidy 4`, the short form of `--sample-ploidy 4`.

```bash
gatk --java-options -Xmx20G HaplotypeCaller \
  -R "$ref" -I "$sample.markdup.bam" \
  -ERC GVCF --sample-ploidy 4 -O "$sample.g.vcf.gz"
```

The archived whole-genome run imports indexed per-accession GVCFs into a GenomicsDB workspace per chromosome. It then runs `GenotypeGVCFs` for each workspace and combines the eight VCFs with `MergeVcfs`. The Methods describe a `CombineGVCFs` consolidation route.

Use `gvcf.list` with one indexed GVCF path per line. The import jobs used `-Xmx200g -Xms200g` and `--bypass-feature-reader true`. The archived joint-genotyping jobs for `chr4.1` and `chr5.1` used `-Xmx450g -Xms450g`. Set the Java memory for the compute allocation and use the same cohort for all eight complete intervals.

```bash
import_heap=${IMPORT_HEAP:-200g}
joint_heap=${JOINT_HEAP:-450g}
inputs=()
for n in {1..8}; do
  chrom="chr${n}.1"
  workspace="merge_${chrom}_database"
  mkdir -p "tmp_$chrom" "${chrom}_tmp"
  gatk --java-options "-Xmx$import_heap -Xms$import_heap" GenomicsDBImport \
    -R "$ref" --variant gvcf.list --intervals "$chrom" \
    --genomicsdb-workspace-path "$workspace" --tmp-dir "tmp_$chrom" \
    --bypass-feature-reader true
  gatk --java-options "-Xmx$joint_heap -Xms$joint_heap" GenotypeGVCFs -R "$ref" -L "$chrom" \
    -V "gendb://$workspace" --sample-ploidy 4 \
    --tmp-dir "${chrom}_tmp" -O "$chrom.cohort.vcf.gz"
  inputs+=(-I "$chrom.cohort.vcf.gz")
done
gatk MergeVcfs "${inputs[@]}" -O all.cohort.vcf.gz
```

## 3. Filter variant sites and prepare the GWAS input

The archived final site selection uses `QUAL > 50`, a positive-depth call count divided by 385 greater than 0.90 (at least 347 positive-depth calls), and more than ten nonreference carriers. The manuscript describes alternate-allele support in at least ten individuals. The counter skips missing sample fields and fields with missing or zero `DP`, then counts a carrier for each remaining `GT` containing a nonreference allele. At a biallelic site, ALT dosage 1, 2, 3 or 4 contributes one carrier per accession. Preserve missing genotypes and the four-allele `GT` field through filtering.

The final predicate operates on the site statistics table, where columns 3–5 contain site quality, positive-depth call count and nonreference carrier count. Use the retained chromosome and position keys to subset the cohort VCF.

```bash
awk 'NR > 1 && $3 > 50 && $4 / 385 > 0.9 && $5 > 10 { print $1 "\t" $2 }' \
  all.merge.vcf.f.gz.info > retained-sites.tsv
```

The final QC call set contains 48,542,275 variants, including 37,112,520 SNPs and 11,429,755 indels. For association testing, match the 358 accessions to the phenotype panel and prepare the SNP dosage matrices as described in [Step 12](12-gwas.md). Record the accession IDs and per-chromosome record counts before and after site filtering.

## Repository example

[`Demos/raw.vcf.gz`](../Demos/raw.vcf.gz) is an example VCF from a subregion of `chr7.1`. [`03.filter_vcf.py`](../03.filter_vcf.py) illustrates filtering by positive depth, site quality and nonreference carrier count. Its fixed carrier cutoff is greater than four. The script expects unphased `GT:AD:DP` fields in that order and emits a simplified GT-only example and a per-site count table.

```bash
python3 03.filter_vcf.py Demos/raw.vcf.gz 0.1 > demo.filtered.vcf
```

## References

- [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360036889312-HaplotypeCaller)
- [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360036727931-GenotypeGVCFs)
