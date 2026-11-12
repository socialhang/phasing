# Step 13 Candidate gene prioritization

Candidate genes were nominated from the trait-associated regions identified by dosage-resolved GWAS in 358 autotetraploid alfalfa accessions. Supplementary Table 43 reports the associated regions and candidate-gene records for the 28 agronomic traits.

## 1. Identify supported association regions

Use the association regions from [Step 12](12-gwas.md). The Methods specify a genome-wide threshold of `P < 1 × 10⁻⁷` (`−log₁₀(P) > 7`) and support from at least three significant neighboring SNPs in the local LD block. Use the reported region coordinates to summarize each association.

## 2. Nominate and annotate candidate genes

Annotate variant effects with SnpEff v5.2a on the same reference used for association testing. The archived `chr.dantao` database uses gene models transferred to the calling reference with Liftoff v1.6.3. Use the gene ID cross-reference to report the corresponding XJDYt2t gene IDs. Coordinate conversion between assemblies requires a separate mapping.

For each retained region, prioritize genes with nonsynonymous protein-coding variants and record the variant effect. The archived candidate scripts rank annotated effects and use the minimum P value across models for each trait and variant. Track this candidate-variant P value separately from the lead SNP P value for the region. Retain separate records when multiple candidates are reported for one region or a gene is reported for multiple traits.

Cross-reference the candidates with *Arabidopsis thaliana* homologues to interpret known gene functions in the context of the associated trait. The reference gene set was functionally annotated using InterProScan v5.77-108.0 against Swiss-Prot, GO, KEGG and InterPro. For example, homologues of the cold-response regulators CBF2 and CBF3 co-localized with winter injury associations.

## 3. Report the candidate set

Supplementary Table 43 records the phenotype, chromosome, region bounds, lead-variant position and P value, candidate-variant position and alleles, variant effect, and gene identifier and name where available. The table contains 576 distinct trait-associated regions and 326 candidate-gene records, representing 298 distinct gene identifiers. Figure 4a shows the association results, and Supplementary Fig. 18 provides supporting plots.

Use the table as the record of the reported candidate set when comparing a reanalysis. The variant-position column identifies candidate variants for follow-up.

## 4. Select candidates for functional follow-up

The study selected six candidates carrying nonsynonymous coding variants from the trait categories chosen for follow-up. *MsLS1* was associated with stem-to-leaf ratio, *MsLS2* with leaf size, and *MsCCX2*, *MsTINY2*, *MsST1* and *MsST2* with salt tolerance. Overexpression assays evaluated projected leaf area for *MsLS1* and *MsLS2* and root elongation under salt stress for the four salt-tolerance candidates (Fig. 4b–d and Extended Data Fig. 4).
