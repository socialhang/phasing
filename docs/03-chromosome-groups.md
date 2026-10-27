# Step 3 Manual partitioning into chromosome groups

The assembly graph was visualized in **Bandage v0.9.0** to assess sequence continuity and manually partition the assembly into **eight homologous chromosome groups**. This was a manual operation; no automated partitioning script was used.

Each group represents one chromosome type and its four homologous copies in the autotetraploid genome. The eight groups are the units for subsequent within-group phasing using F1 progeny data. This grouping step does not itself resolve the four haplotypes within each group.

## Reference alignment supporting manual grouping

The recorded alignment command is:

```bash
minimap2 -x asm5 -t 40 fll5gp.chr.p.60.fa hifiasm.asm.bp.p_utg.fa > scm.paf
```

| Item | Recorded value |
| --- | --- |
| Reference/target FASTA | `fll5gp.chr.p.60.fa` |
| Query unitig FASTA | `hifiasm.asm.bp.p_utg.fa` |
| Aligner | minimap2; historical version not yet confirmed |
| Preset | `-x asm5` |
| Threads | `-t 40` |
| Output | `scm.paf` |
| Original PBS allocation | `high1` queue, one node, 40 CPUs, `mem=76g` |

The reference filename and its link to an earlier project assembly were verified. Its formal assembly name, release version, and public accession have not yet been established; the filename should not be treated as an accession. The unitig FASTA was linked from the `02.hifiasm_ont` assembly directory.

The [portable alignment script](../scripts/03_reference_alignment.sh) preserves the original preset, thread count, and reference/query order. For example, include the following in an appropriate PBS job:

```bash
bash /path/to/phasing/scripts/03_reference_alignment.sh \
  /path/to/fll5gp.chr.p.60.fa \
  /path/to/hifiasm.asm.bp.p_utg.fa \
  /path/to/scm.paf
```

`MINIMAP2` can specify the executable path. This command produces supporting alignments; the graph partition itself remains manual.

### Assignment records found in the analysis directory

| Record | Observed contents or role |
| --- | --- |
| `scm.paf` | Alignment output specified by the original command. |
| `scm.best.paf` | A selected alignment file; its generation command and selection criteria were not found in the inspected scripts. |
| `utg2chr.txt` | Two whitespace-separated columns: unitig ID and chromosome-group label. |
| `chr1.utgid` to `chr8.utgid` | Per-group node ID lists. |
| `chr1.a.utgid` to `chr8.a.utgid` | A second set of per-group node lists; the suffix alone does not establish which lists were used for the final graph partition. |
| `utg2chr.addtelo.txt` | A comma-separated annotation file with header `utgid,chrid`. |

For example, an inspected alignment for `utg000001l` targets reference sequence `chr1.2`, while the corresponding row in `utg2chr.txt` assigns it to `chr1`. This documents an example of a reference sequence assignment being represented at chromosome-group level; it does not establish the general tie-breaking or label-conversion rule. The assignment table also includes a `Chloroplast` label, which is distinct from the eight nuclear chromosome groups.

The supporting [add_chr.py helper](../scripts/chromosome_groups/add_chr.py) prints the nodes in a supplied list, followed by nodes assigned to the matching group in `utg2chr.txt` that are absent from the supplied list. It does not infer assignments from PAF, edit the graph, or perform haplotype phasing. Its original logic is retained. To use it, work in the directory containing `utg2chr.txt` and supply a basename such as `chr1.utgid`, since the helper derives the group name from the text before the first period in the argument:

```bash
# Example use of the recovered helper; the output filename is illustrative.
python3 /path/to/phasing/scripts/chromosome_groups/add_chr.py chr1.utgid \
  > chr1.supplemented.utgid
```

## Inputs and manual workflow

1. Align the unitig FASTA to the reference using the recorded minimap2 command and inspect the chromosome-assignment evidence.
2. Open the assembly GFA in Bandage v0.9.0.
3. Inspect the graph representation and assess sequence continuity.
4. Use the reference-alignment chromosome assignments for most nodes as the primary evidence for manually partitioning the graph into eight homologous chromosome groups.
5. Use the resulting group assignments as the basis for subsequent F1 progeny-based phasing within each chromosome group.

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

Record unresolved or shared nodes explicitly rather than assigning them by an undocumented default. Screenshots of each group and any exported subgraphs can accompany the table. The original directory contains assignment tables and node lists, as described above; their final selection and correspondence to the manually partitioned graphs still need to be recorded. The three-column format above is a recommended audit record, not a description of the original two-column table.

The following details remain to be added from the original manual procedure:

- The exact GFA used and its relationship to the graph evaluated in Step 2.
- The formal identity and release/accession of the reference FASTA.
- The historical minimap2 version, the selection rule for `scm.best.paf`, and the criteria for assigning alignments to chromosome-group labels.
- How disconnected fragments, repeat-associated connections, and ambiguous/shared nodes were handled.
- Whether edges or nodes were edited, removed, or duplicated during grouping, and how those actions were recorded.
- The format of the assignments or group-specific graph files passed to the next step.

## Methods text

> Unitig sequences were aligned to the reference assembly (`fll5gp.chr.p.60.fa`) using minimap2 with the parameters `-x asm5 -t 40`. The assembly graph was subsequently visualized using Bandage (v0.9.0) to assess sequence continuity and manually partition the assembly into eight homologous chromosome groups, guided primarily by the reference-alignment chromosome assignments for most nodes. Each group comprised the four homologous copies of one chromosome type and was used as a separate unit for subsequent phasing with F1 progeny data. Under random chromosome segregation without double reduction or preferential pairing, the joint probability that two specified homologous chromosome copies are transmitted to the same gamete is 1/6, whereas the corresponding probability for specified copies from independently assorting chromosome groups is 1/4. Assigning graph fragments to their respective chromosome groups therefore provides a framework for interpreting progeny segregation patterns during subsequent phasing.

This paragraph records the supplied description with explicit model assumptions. The reference genome, alignment details, and handling of unresolved nodes listed above are needed to complete the reproducibility description.

## References

- [minimap2 source and documentation](https://github.com/lh3/minimap2)
- [Bandage source and documentation](https://github.com/rrwick/Bandage)
- [Bandage releases](https://github.com/rrwick/Bandage/releases)
- [Theoretical basis for genetic linkage analysis in autotetraploid species](https://pmc.ncbi.nlm.nih.gov/articles/PMC406462/)
