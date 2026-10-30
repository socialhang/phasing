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
| Aligner | minimap2 |
| Preset | `-x asm5` |
| Threads | `-t 40` |
| Output | `scm.paf` |
| Original PBS allocation | `high1` queue, one node, 40 CPUs, `mem=76g` |

The reference FASTA is `fll5gp.chr.p.60.fa`. The query unitig FASTA is from the `02.hifiasm_ont` assembly directory.

The [portable alignment script](../scripts/03_reference_alignment.sh) preserves the original preset, thread count, and reference/query order. For example, include the following in an appropriate PBS job:

```bash
bash /path/to/phasing/scripts/03_reference_alignment.sh \
  /path/to/fll5gp.chr.p.60.fa \
  /path/to/hifiasm.asm.bp.p_utg.fa \
  /path/to/scm.paf
```

`MINIMAP2` can specify the executable path. This command produces supporting alignments; the graph partition itself remains manual.

### Chromosome-assignment files

| File | Contents or role |
| --- | --- |
| `scm.paf` | Alignment output specified by the original command. |
| `scm.best.paf` | Selected alignments. |
| `utg2chr.txt` | Two whitespace-separated columns: unitig ID and chromosome-group label. |
| `chr1.utgid` to `chr8.utgid` | Per-group node ID lists. |
| `chr1.a.utgid` to `chr8.a.utgid` | Additional per-group node lists. |
| `utg2chr.addtelo.txt` | A comma-separated annotation file with header `utgid,chrid`. |

For example, an alignment for `utg000001l` targets reference sequence `chr1.2`, while the corresponding row in `utg2chr.txt` assigns it to `chr1`. This represents the reference sequence assignment at chromosome-group level. The assignment table also includes a `Chloroplast` label, which is distinct from the eight nuclear chromosome groups.

The supporting [add_chr.py helper](../scripts/chromosome_groups/add_chr.py) prints the nodes in a supplied list, followed by nodes assigned to the matching group in `utg2chr.txt` that are absent from the supplied list. To use it, work in the directory containing `utg2chr.txt` and supply a basename such as `chr1.utgid`, since the helper derives the group name from the text before the first period in the argument:

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

The grouping was guided primarily by reference-genome alignment assignments for most nodes. The reference provides chromosome-group assignment evidence; this step does not assign the four homologous copies to separate haplotypes.

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

## Methods text

> Unitig sequences were aligned to the reference assembly (`fll5gp.chr.p.60.fa`) using minimap2 with the parameters `-x asm5 -t 40`. The assembly graph was subsequently visualized using Bandage (v0.9.0) to assess sequence continuity and manually partition the assembly into eight homologous chromosome groups, guided primarily by the reference-alignment chromosome assignments for most nodes. Each group comprised the four homologous copies of one chromosome type and was used as a separate unit for subsequent phasing with F1 progeny data. Under random chromosome segregation without double reduction or preferential pairing, the joint probability that two specified homologous chromosome copies are transmitted to the same gamete is 1/6, whereas the corresponding probability for specified copies from independently assorting chromosome groups is 1/4. Assigning graph fragments to their respective chromosome groups therefore provides a framework for interpreting progeny segregation patterns during subsequent phasing.

## References

- [minimap2 source and documentation](https://github.com/lh3/minimap2)
- [Bandage source and documentation](https://github.com/rrwick/Bandage)
- [Bandage releases](https://github.com/rrwick/Bandage/releases)
- [Theoretical basis for genetic linkage analysis in autotetraploid species](https://pmc.ncbi.nlm.nih.gov/articles/PMC406462/)
