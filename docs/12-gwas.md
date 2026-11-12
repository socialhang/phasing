# Step 12 Dosage-resolved genome-wide association analysis

The diversity-panel analysis assessed 28 agronomic traits in 358 autotetraploid alfalfa accessions. The dosage scan covers all eight complete chromosomes of the H1 calling reference, `chr1.1` through `chr8.1`. The study reports 22,256,624 variant sites for the GWAS. The short `chr7.1` VCF in `Demos/` illustrates variant processing. The 100-progeny k-mer association analysis for *ML1* on chromosome 5 is a separate analysis.

## Inputs

| Input | Contents |
| --- | --- |
| Filtered genotypes | Tetraploid genotype calls from [Step 11](11-variant-calling-qc.md), converted to a numeric dosage matrix for association testing. |
| Reference and coordinate dictionary | The eight chromosomes of the reference used for variant calling. Keep that reference and its chromosome names with the genotype matrix, score files and reported intervals. |
| Sample and trait records | Matched genotype and phenotype identifiers, a trait-name crosswalk, and accession-level values for each trait. The phenotype merge represents unavailable records as `NA`. |

For biallelic VCF records, encode ALT dosage as the number of allele `1` entries in a complete tetraploid `GT`, from 0 to 4. Spot checks of the archived numeric matrix against the cohort VCF confirm ALT counting. Store missing calls as `NA`. The original genotype-export helper skipped sites with any allele longer than one base and recoded calls with `*` or missing alleles to `NA`. Match the 358 sample IDs to the phenotype file and record the number of accessions with usable values per trait. Preserve the variant's chromosome, position and allele orientation when preparing the dosage tables.

## Association model and chromosome coverage

The archived numeric drivers import one trait and its genotypes into GWASpoly with `ploidy = 4`, `format = "numeric"` and `n.traits = 1`. They request a leave-one-chromosome-out (LOCO) kinship calculation with `set.K(..., LOCO = TRUE)` and use the maximum genotype-frequency setting `set.params(geno.freq = 1 - 5/385)`. The retained numeric driver supplies no principal components or other fixed effects.

The archived job scripts dispatch a trait across `chr1.1.csv` through `chr8.1.csv` and collect chromosome-specific results. For a reproducible LOCO analysis, preserve the genome-wide marker input used to construct each kinship matrix and record which chromosome was excluded. The archived numeric driver requests LOCO while the job scripts pass chromosome-specific genotype files. These files establish the requested option and analysis coverage, but do not establish the effective kinship matrices used for the retained scores.

The manuscript lists six tests. The archived score tables used for downstream candidate selection contain those six plus two duplex-dominance columns:

| Tests | Score columns |
| --- | --- |
| Additive and general | `additive`, `general` |
| Diploidized | `diplo-additive`, `diplo-general` |
| Simplex dominance | `1-dom-alt`, `1-dom-ref` |
| Duplex dominance, present in the archived scores | `2-dom-alt`, `2-dom-ref` |

The retained file named `run.single_allmodel.R` calls an additive-only scan. The eight-column score tables are the starting point for reproducing the downstream catalogue. Reproducing their association fits requires the matching full-model driver, package version and kinship inputs.

For a new full-chromosome LOCO scan, use a single dosage table with markers from all eight chromosomes so that each kinship calculation can exclude the chromosome under test. Give the table three leading columns for marker ID, chromosome and position, then one dosage column per accession. The phenotype CSV has accession ID first and one named trait column. This example requests all model columns present in the archived score tables.

```r
library(GWASpoly)

threads <- as.integer(Sys.getenv("GWAS_THREADS", "1"))
data <- read.GWASpoly(
  ploidy = 4, pheno.file = "trait.csv", geno.file = "all-chromosomes.dosage.csv",
  format = "numeric", n.traits = 1, delim = ","
)
kinship <- set.K(data, LOCO = TRUE, n.core = threads)
params <- set.params(geno.freq = 1 - 5/385)
models <- c("additive", "diplo-general", "diplo-additive", "general", "1-dom", "2-dom")
fit <- GWASpoly(kinship, models = models, params = params, n.core = threads)
scored <- set.threshold(fit, method = "Bonferroni", level = 0.05)
trait <- names(read.csv("trait.csv", nrows = 0, check.names = FALSE))[2]
write.GWASpoly(scored, trait = trait, filename = "association-scores.tsv", what = "scores", delim = "\t")
```

The Bonferroni setting reproduces the driver's plotting threshold. Apply the locus-reporting criterion below to the exported association P values.

## Combine scores across chromosomes

GWASpoly exports marker-level scores by chromosome. The archived score-merging helper joins the chromosome files for each trait, converts each score `s` to a P value as `10^(-s)`, and preserves the eight model columns. For this export, a missing score is represented by `s = 0`, and a row is omitted if its highest score is below 2. Keep the unfiltered scores for Q–Q diagnostics and for the REF and ALT fields, which are absent from the merged export. Use the same trait IDs across files.

## Report trait-associated loci

The manuscript uses a lead-SNP threshold of `P < 1 × 10⁻⁷` and describes support by at least three neighboring SNPs with `P < 1 × 10⁻⁷` within an approximately 8-kb local LD block. The threshold was motivated by the reported LD decay distance and an effective-test estimate. Treat it as the manuscript’s reporting rule.

The retained catalogue contains 576 distinct trait–chromosome–interval records across 28 traits. The downstream candidate workflow uses reviewed trait-specific intervals and, for candidate lookup, the minimum P value among the columns in the retained score table. A historical locus-finding helper merged successive signals at `P < 1 × 10⁻⁶` separated by at most 100 kb and retained groups containing a lead signal at `P < 1 × 10⁻⁷`. That helper documents an earlier exploratory rule. Use the retained interval table to reproduce the reported loci.

Inspect unfiltered Q–Q plots by trait and model. Residual inflation and the imperfect correspondence between historical accession-level phenotypes and genotypes from sequenced plants affect interpretation of these associations. Retain the association model, sample count, lead marker and its P value when reporting a locus. See [Step 13](13-candidate-prioritization.md) for functional annotation and candidate selection within the intervals.

## Outputs

| Output | Contents |
| --- | --- |
| Chromosome-wise scores | Marker identifiers, positions and model-specific scores for all eight H1 chromosomes. |
| Trait-wise association tables | Combined model-specific P values with the export filtering recorded. |
| Q–Q and Manhattan plots | Trait- and model-specific diagnostics and genomic signal distributions. |
| Trait-associated intervals | Reviewed locus coordinates and lead markers for candidate prioritization. |

## Reference

- [Rosyara et al., Software for Genome-Wide Association Studies in Autopolyploids and Its Application to Potato](https://doi.org/10.3835/plantgenome2015.08.0073)
- [GWASpoly package vignette](https://jendelman.github.io/GWASpoly/GWASpoly.html)
