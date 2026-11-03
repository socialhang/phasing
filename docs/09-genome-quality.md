# Step 9 Genome quality assessment

The final XJDYt2t and Adina assemblies were assessed for sequence continuity, read support, consensus accuracy, gene-space completeness, and haplotype consistency. Progeny inheritance signals and Pore-C contacts were evaluated separately for switch errors. Hi-C and Pore-C contact maps were inspected for large-scale structural misassemblies.

The [evaluation wrapper](../scripts/09_genome_quality.sh) provides independent stages. Run each stage in its own output directory and submit computational jobs through PBS. Examples use `assembly.fa` for the complete assembly and `hap1.fa` for one haplotype. Nuclear chromosome identifiers follow `chr1.1` through `chr8.4`: the first number denotes the chromosome group and the second the haplotype. Apply the same procedure to both parental assemblies and all four haplotypes.

## Software and settings

| Component | Version and settings |
| --- | --- |
| [BWA-MEM](https://github.com/lh3/bwa) | 0.7.17-r1198; default alignment settings, `-t 40` |
| [Winnowmap](https://github.com/marbl/Winnowmap) | 2.03; `-ax map-pb` for HiFi, `-ax map-ont` for UL reads |
| [Minimap2](https://github.com/lh3/minimap2) | `-ax map-hifi` for HiFi, `-ax map-ont` for UL reads; 40 threads |
| [GCI](https://github.com/yeeus/GCI) | 1.0; `-mq 0`, 10 threads |
| [CRAQ](https://github.com/JiaoLaboratory/CRAQ) | 1.10; `-q 0`, 10 threads; other analysis settings at defaults |
| [Merqury](https://github.com/marbl/merqury) | 1.3; NGS k-mer database built with `k=21` |
| [BUSCO](https://busco.ezlab.org/busco_userguide.html) | 5.4.7; `embryophyta_odb10`, 10 threads |
| [HiGlass](https://docs.higlass.io/) | 1.13.4; structural inspection of Hi-C and Pore-C contact maps |

Additional command-line dependencies are samtools, bedtools, meryl, paftools.js, pigz, Python 3, Perl, and mosdepth for optional depth plots. The processing helpers use the Python standard library; the optional support-separation plotting helper requires NumPy and Matplotlib.

## 1 Basic assembly statistics

Calculate total sequence length, sequence N50, contig N50, and gap counts:

```bash
bash /path/to/phasing/scripts/09_genome_quality.sh basic assembly.fa > assembly.metrics.json
```

The [statistics helper](../scripts/quality/assembly_metrics.py) defines a gap as a maximal run of one or more `N` bases. Contigs are the nonempty sequences between these runs. Contig N50 is the length at which the cumulative length of contigs sorted from longest to shortest reaches at least half of their total length. Sequence N50 is also reported, before splitting at gaps.

Inspect both chromosome ends for terminal plant telomeric repeats (`CCCTAAA` and its reverse complement `TTTAGGG`). Record the reviewed ends in a tab-separated file, with one row per chromosome end:

```text
# chromosome  end    telomere_present
chr1.1        left   1
chr1.1        right  1
```

Use `0` for an end without an identified telomere and `1` for an end with a telomere. Supply all reviewed ends to count telomere-positive chromosome ends:

```bash
bash /path/to/phasing/scripts/09_genome_quality.sh basic assembly.fa \
  --telomere-ends telomere_ends.tsv > assembly.metrics.json
```

## 2 Short-read mapping and coverage

Align paired-end WGS reads to the complete assembly with BWA-MEM, sort and index the alignments, and summarize mapping rates and per-sequence coverage:

```bash
bwa index assembly.fa
bwa mem -t 40 assembly.fa ngs_R1.fq.gz ngs_R2.fq.gz \
  | samtools sort -@ 10 -o ngs.bam
samtools index ngs.bam
samtools flagstat ngs.bam > ngs.flagstat.txt
samtools coverage ngs.bam > ngs.coverage.tsv
```

The wrapper's `ngs` stage runs these commands. `flagstat` reports alignment and pairing statistics; `coverage` reports covered bases, coverage fraction, and mean depth for each sequence.

## 3 Haplotype-specific long-read support

First align HiFi and UL reads independently to the complete assembly with Winnowmap. Prepare the repetitive 15-mer list as described in the Winnowmap workflow:

```bash
meryl count k=15 output repetitive.meryl assembly.fa
meryl print greater-than distinct=0.9998 repetitive.meryl > repetitive_k15.txt
winnowmap -W repetitive_k15.txt -t 40 -ax map-pb assembly.fa hifi.fq.gz > hifi.w.sam
winnowmap -W repetitive_k15.txt -t 40 -ax map-ont assembly.fa ul.fq.gz > ul.w.sam
paftools.js sam2paf hifi.w.sam > hifi.w.paf
paftools.js sam2paf ul.w.sam > ul.w.paf
```

The wrapper's `winnowmap` stage additionally creates sorted, indexed BAM files. Its `partition` stage retains the first mapper-ranked alignment for each read, matching the source `best0 ... 1` selection, and assigns the read according to the target chromosome's haplotype suffix. Keep PAF records in mapper order for this selection. Reads whose selected target is outside the four nuclear haplotypes are excluded.

Create `hifi_reads.list` and `ul_reads.list`, each containing one gzip-compressed FASTQ path per line. In separate directories for HiFi and UL, run:

```bash
bash /path/to/phasing/scripts/09_genome_quality.sh partition /path/to/hifi.w.paf hifi_reads.list hifi
bash /path/to/phasing/scripts/09_genome_quality.sh partition /path/to/ul.w.paf ul_reads.list ul
```

Each run produces `hap1.reads.txt` through `hap4.reads.txt` and the corresponding extracted FASTQs. Remap each read set to its respective haplotype using Minimap2:

```bash
bash /path/to/phasing/scripts/09_genome_quality.sh remap hap1.fa hap1.hifi.fq.gz hifi hifi.hap1
bash /path/to/phasing/scripts/09_genome_quality.sh remap hap1.fa hap1.ul.fq.gz ul ul.hap1
```

These commands use `-ax map-hifi` and `-ax map-ont`, respectively, and generate SAM, sorted/indexed BAM, and target-coordinate-sorted PAF files. Repeat for haplotypes 2–4.

## 4 GCI and effective contig N50

GCI evaluates assembly continuity using the haplotype-specific HiFi and UL alignments:

```bash
GCI -r hap1.fa \
  --hifi hifi.hap1.bam hifi.hap1.paf \
  --nano ul.hap1.bam ul.hap1.paf \
  --chrs chr1.1,chr2.1,chr3.1,chr4.1,chr5.1,chr6.1,chr7.1,chr8.1 \
  -t 10 -p -it pdf -mq 0 -d gci.hap1
```

The corresponding wrapper stage is `gci hap1.fa 1 hifi.hap1.bam hifi.hap1.paf ul.hap1.bam ul.hap1.paf`. Retain the GCI scores and uncovered-region BED output. To calculate effective contig N50, mask those uncovered intervals with `N` bases and calculate contig N50 after splitting at the resulting runs of `N`:

```bash
bash /path/to/phasing/scripts/09_genome_quality.sh effective hap1.fa uncovered_regions.bed
```

The `contig_N50_bp` field in `effective.metrics.json` is the effective contig N50. Use the uncovered-region BED for the evaluated haplotype from the GCI run; BED coordinates are zero-based, half-open.

## 5 CRAQ assembly error assessment

Analyze the short-read and long-read BAM alignments with CRAQ to estimate assembly errors and calculate regional and structural assembly quality indices, R-AQI and S-AQI:

```bash
craq -ngs ngs.bam -sms ul.hap1.bam -g hap1.fa -q 0 -t 10 -D craq.hap1
```

The wrapper stage is `craq hap1.fa ngs.bam ul.hap1.bam craq.hap1`. BAM reference identifiers and lengths must match the evaluated sequences. Repeat for each haplotype and retain both CRAQ summary statistics and error intervals.

For optional 100-kb depth profiles, the `depth` stage runs:

```bash
mosdepth -n --fast-mode --by 100000 -t 40 hap1.depth ul.hap1.bam
```

## 6 Merqury completeness and consensus QV

Build a 21-mer database from the short-read WGS data and evaluate the complete assembly:

```bash
meryl count k=21 output ngs.k21.meryl ngs_R1.fq.gz ngs_R2.fq.gz
bash "$MERQURY/merqury.sh" ngs.k21.meryl assembly.fa assembly
```

Set `MERQURY` to the Merqury v1.3 installation directory. The wrapper stage is `merqury assembly.fa ngs_R1.fq.gz ngs_R2.fq.gz`. The k-mer length is set during database construction, rather than passed as a `-k` argument to `merqury.sh`. Summarize Merqury's read-k-mer completeness, consensus QV, and estimated base error rate. For the tetraploid assembly, `assembly.fa` contains all four haplotypes.

## 7 BUSCO gene-space completeness

Assess each haplotype separately using BUSCO v5.4.7 and `embryophyta_odb10`:

```bash
busco -m genome -i hap1.fa -o busco.hap1 \
  -l /path/to/embryophyta_odb10 -c 10 --offline
```

Set `BUSCO_LINEAGE` to the local lineage directory when using the wrapper's `busco genome hap1.fa busco.hap1` stage. For evaluation of predicted protein sets, use `busco proteins hap1.pep busco.hap1.proteins`. Report complete single-copy, complete duplicated, fragmented, and missing BUSCOs, distinguishing genome-mode results from protein-mode results.

## 8 Switch errors from progeny and Pore-C data

### Progeny correlation support

Identify genome-wide unique 51-mers in the final assembly and retain one randomly selected marker per informative 10-kb window. Construct its presence/absence vector across the 100 progeny using the frequency threshold **at least three** described in [Step 4](04-progeny-phasing.md). Approximately 310,000 informative windows give approximately 96 billion ordered window-pair comparisons, including the diagonal before filtering.

Calculate pairwise Pearson correlations between the vectors and encode the correlation signals as the pseudo-contact records described in Step 4. The matrix counts represent correlation-derived support. For each focal window, accumulate that support separately for the four haplotypes of the same chromosome group, then identify the haplotype with the greatest accumulated support. Thus, the comparison uses the strongest **haplotype-level correlation support**.

### Shared switch-region filtering

The public workflow follows the directory scripts:

1. Convert contact positions to 10-kb window indices using integer division by 10,000; discard pairs within the same window.
2. Count repeated ordered window pairs to construct the window matrix.
3. Accumulate matrix support by the partner window's haplotype, considering only pairs within the same chromosome group. Both endpoints contribute to their respective window summaries.
4. Retain windows whose maximum-support haplotype differs from their assembled haplotype and whose summed support across all four haplotypes is **greater than 200**.
5. Convert candidate windows to 10-kb BED intervals. Merge intervals separated by **at most 50,000 bp**, then retain merged regions of length **greater than 10,000 bp**.
6. Divide the total retained merged-region length by the full length of the evaluated haplotype. Report progeny-derived and Pore-C-derived rates separately.

For Pore-C, replace the progeny pseudo-contact matrix with the Pore-C contact matrix and apply the same accumulation and interval filtering. Supply MND-format records for either modality, with chromosome and position in columns 2–3 and 6–7:

```text
strand1 chr1 pos1 fragment1 strand2 chr2 pos2 fragment2 ...
```

Create the assembly FASTA index, then run the two modalities independently:

```bash
samtools faidx assembly.fa
bash /path/to/phasing/scripts/09_genome_quality.sh switch progeny.mnd assembly.fa.fai progeny
bash /path/to/phasing/scripts/09_genome_quality.sh switch porec.mnd assembly.fa.fai porec
```

Outputs include window-pair counts (`*.cor.tsv`), per-window haplotype support (`*.stat.tsv`), candidate windows, retained merged switch regions, and per-haplotype/overall switch error rates. The statistics columns are: window ID, assembled haplotype, maximum-support haplotype, and support for haplotypes 1–4. Merged-region lengths include the intervening spans joined by the 50-kb merge rule; isolated 10-kb candidates are removed by the length filter.

The [window converter](../scripts/quality/window_contacts.py) and [support accumulator](../scripts/quality/sum_haplotype_support.py) retain the source implementations. Their detailed behavior is described in the [helper notes](../scripts/quality/README.md).

### Visualize strongest and second-strongest support

The [plotting helper](../scripts/quality/plot_support_separation.py) reads the seven-column `*.stat.tsv` tables before candidate selection and interval merging. For each window, let `T` be the sum of the four haplotype supports and `S1 >= S2` their two highest values. It plots a hexbin density map of `S1/T` versus `S2/T` and per-haplotype empirical cumulative distributions of `(S1-S2)/T`, including a low-margin inset. These quantities are normalized accumulated supports, not individual Pearson coefficients.

```bash
python3 /path/to/phasing/scripts/quality/plot_support_separation.py \
  --dataset 'XJDYt2t progeny=progeny.stat.tsv' \
  --dataset 'XJDYt2t Pore-C=porec.stat.tsv' \
  --outdir support_separation
```

By default, only windows with `T > 200` enter the plots. Outputs include PNG/PDF/SVG figures, per-window metrics, per-haplotype summary counts, and input SHA-256 checksums. Exactly tied windows have zero normalized margin. The summary also reports fractions with margins at most 0.05 and 0.1 as descriptive statistics; these values are not additional switch-error thresholds. Windows absent from the input support tables are not counted as recorded or evaluated windows.

## 9 Structural inspection

Visualize Hi-C and Pore-C contact maps in HiGlass v1.13.4. Inspect each chromosome and compare the four homologous copies for continuity of the diagonal and anomalous off-diagonal blocks. Review candidate switch regions alongside these maps and the progeny support summaries. Record the inspected assembly version, contact-map resolution, and any reviewed intervals with the final assessment outputs.

## Methods text

Basic assembly statistics included contig N50, telomere counts, and gap counts. Short-read WGS data were aligned using BWA-MEM v0.7.17-r1198 with default alignment parameters to evaluate mapping rates and coverage. HiFi and UL reads were aligned using Winnowmap v2.03 (`-x map-pb` and `-x map-ont`, respectively), assigned to haplotypes according to their selected best alignments, and remapped to the corresponding haplotypes using Minimap2 (`-x map-hifi` and `-x map-ont`). GCI v1.0 (`-mq 0`) was used to calculate assembly continuity scores and identify uncovered regions; effective contig N50 was calculated after splitting the assembly at these regions. CRAQ v1.10 (`-q 0`) was used to estimate assembly errors and calculate R-AQI and S-AQI. Merqury v1.3, using a 21-mer database constructed from NGS reads, assessed k-mer completeness, consensus QV, and base error rates. Gene-space completeness was evaluated separately for each haplotype using BUSCO v5.4.7 with `embryophyta_odb10`.

For switch-error evaluation, one unique 51-mer was randomly selected per informative 10-kb window, and its segregation vector across 100 progeny was used to calculate pairwise Pearson correlations. Correlation-derived support in the pseudo-contact matrix was accumulated for each of the four homologous haplotypes. Windows were identified as candidates when their maximum-support haplotype differed from their assembled haplotype and their total support exceeded 200. Candidate intervals separated by at most 50 kb were merged, and merged regions of length at most 10 kb were removed. Switch error rates were calculated as the retained merged-region length divided by the full haplotype length. Pore-C switch error rates were evaluated with the same procedure using the Pore-C contact matrix. Hi-C and Pore-C contact maps were visualized using HiGlass v1.13.4 for structural inspection.
