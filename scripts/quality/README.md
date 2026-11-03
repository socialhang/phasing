# Quality-assessment helpers

The entry point is [09_genome_quality.sh](../09_genome_quality.sh); see the [Step 9 guide](../../docs/09-genome-quality.md) for commands and interpretation.

## Source implementations

- `extract_reads_from_list.py`: FASTQ extraction from a list of gzip-compressed inputs; retains the original four-line FASTQ processing and read-name matching.
- `window_contacts.py`: maps both endpoints of MND records to 10-kb windows and excludes identical-window pairs.
- `sum_haplotype_support.py`: accumulates support for the four homologous haplotypes within each chromosome group and reports the maximum-support haplotype.
- `count_window_pairs.pl`: equivalent to the source `count0` line counter; counts complete, ordered window-pair records.

The source support accumulator initializes a window's four counters on its first encountered pair and adds support on subsequent encounters. The initializing contribution is therefore omitted. Ties select the first maximum, in haplotype order 1–4. The Perl counter emits records in hash order, so that initializing pair can depend on record order. These source behaviors are retained for comparison with the directory workflow. Matrix counts are correlation-derived pseudo-contact support for progeny, or contact support for Pore-C; they are accumulated numerically without converting them back to individual Pearson coefficients.

## Portable helpers

- `assign_haplotypes.py`: first mapper-ranked PAF alignment per read, equivalent to the source `best0 ... 1`, followed by grouping on chromosome suffix `.1`–`.4`. Input must retain mapper order.
- `assembly_metrics.py`: sequence/contig N50, counts and lengths of N-runs, and optional counting of reviewed telomere-positive chromosome ends. It reuses the repository FASTA reader.
- `switch_summary.py`: retained, merged interval length divided by the full nuclear haplotype length from a FASTA index. BED intervals must be nonoverlapping. The overall rate uses the sum of all four haplotype lengths as its denominator.
- `plot_support_separation.py`: optional NumPy/Matplotlib plotting of normalized strongest/second-strongest supports and empirical support-margin distributions, with a low-margin inset. Exports vector figures, window metrics, counts, and input checksums; reads existing support tables without recalculating or modifying them.

Chromosome identifiers must follow `chr1.1` through `chr8.4`. Keep assembly identifiers identical across alignments, contact records, and the FASTA index. Genomic analyses should be submitted through the cluster scheduler.
