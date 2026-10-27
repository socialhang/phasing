# Step 3 Manual partitioning into chromosome groups

The assembly graph was visualized in **Bandage v0.9.0** to assess sequence continuity and manually partition the assembly into **eight homologous chromosome groups**. This was a manual operation; no automated partitioning script was used.

Each group represents one chromosome type and its four homologous copies in the autotetraploid genome. The eight groups are the units for subsequent within-group phasing using F1 progeny data. This grouping step does not itself resolve the four haplotypes within each group.

## Inputs and manual workflow

1. Open the assembly GFA in Bandage v0.9.0.
2. Inspect the graph representation and assess sequence continuity.
3. Use the chromosome assignments obtained by aligning most graph nodes to an existing reference genome as the primary evidence for manually partitioning the graph into eight homologous chromosome groups.
4. Use the resulting group assignments as the basis for subsequent F1 progeny-based phasing within each chromosome group.

The grouping was guided primarily by reference-genome alignment assignments for most nodes. The reference provides chromosome-group assignment evidence; this step does not assign the four homologous copies to separate haplotypes. Treatment of nodes with ambiguous, conflicting, or absent reference alignments has not yet been specified.

These steps describe the supplied workflow at its currently documented level. No numerical grouping threshold or automatic classification rule has been supplied. Visual inspection should not be interpreted as proof that all graph connections are correct.

## Segregation rationale

Under a balanced, random chromosome-segregation model without double reduction or preferential pairing, a diploid gamete receives two of the four homologous copies of each chromosome type. For four copies A1, A2, A3, and A4, the six possible unordered pairs are:

```text
A1+A2  A1+A3  A1+A4  A2+A3  A2+A4  A3+A4
```

These pairs have equal probabilities in this model. Thus, the **joint probability that a specified pair of distinct homologous copies is present in the same gamete** is:

```text
P(A1 and A2) = 1 / choose(4, 2) = 1/6
```

For specified copies from two different chromosome groups, independent assortment gives a marginal inclusion probability of 1/2 for each copy and hence:

```text
P(A1 and B1) = (1/2) × (1/2) = 1/4
```

These are joint transmission probabilities for specified chromosome copies from one parent, not conditional probabilities and not the probability that any pair of chromosomes is transmitted. The distinction between homologous and independently assorting chromosome copies motivates separating the graph into chromosome groups before using progeny segregation information for phasing.

These probabilities are theoretical expectations under the stated assumptions. They are not fixed empirical thresholds for fragment assignment. Double reduction, preferential pairing, segregation distortion, recombination at individual marker loci, and marker specificity can affect observed progeny patterns. The probabilities alone do not establish a fragment's chromosome assignment.

## Records needed for an auditable manual partition

To make the manual assignments traceable, retain the exact input GFA and an assignment table such as:

```text
node_id    chromosome_group    assignment_evidence
```

Record unresolved or shared nodes explicitly rather than assigning them by an undocumented default. Screenshots of each group and any exported subgraphs can accompany the table. This is a recommended record format; no project-specific assignment table or screenshots have yet been supplied with this step.

The following details remain to be added from the original manual procedure:

- The exact GFA used and its relationship to the graph evaluated in Step 2.
- The reference genome name, assembly version/accession, and chromosome identifiers.
- The alignment software, version, parameters, and criteria used to assign nodes to reference chromosomes.
- How disconnected fragments, repeat-associated connections, and ambiguous/shared nodes were handled.
- Whether edges or nodes were edited, removed, or duplicated during grouping, and how those actions were recorded.
- The format of the assignments or group-specific graph files passed to the next step.

## Methods text

> Subsequently, the assembly graph was visualized using Bandage (v0.9.0) to assess sequence continuity and manually partition the assembly into eight homologous chromosome groups, guided primarily by chromosome assignments obtained from alignments of most graph nodes to an existing reference genome. Each group comprised the four homologous copies of one chromosome type and was used as a separate unit for subsequent phasing with F1 progeny data. Under random chromosome segregation without double reduction or preferential pairing, the joint probability that two specified homologous chromosome copies are transmitted to the same gamete is 1/6, whereas the corresponding probability for specified copies from independently assorting chromosome groups is 1/4. Assigning graph fragments to their respective chromosome groups therefore provides a framework for interpreting progeny segregation patterns during subsequent phasing.

This paragraph records the supplied description with explicit model assumptions. The reference genome, alignment details, and handling of unresolved nodes listed above are needed to complete the reproducibility description.

## References

- [Bandage source and documentation](https://github.com/rrwick/Bandage)
- [Bandage releases](https://github.com/rrwick/Bandage/releases)
- [Theoretical basis for genetic linkage analysis in autotetraploid species](https://pmc.ncbi.nlm.nih.gov/articles/PMC406462/)
