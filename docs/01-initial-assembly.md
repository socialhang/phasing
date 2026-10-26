# Step 1 Initial assembly

This step documents the initial hifiasm assembly command, its input files, and its explicit parameters. The accompanying [Bash script](../scripts/01_initial_assembly.sh) preserves the analysis options and input order from the original `work.sh`, with a configurable executable path and working directory.

## Software and assembly mode

- **Assembler:** hifiasm **0.25.0-r726**, verified using the executable's `--version` output.
- **Long-read mode:** the ONT simplex assembly preset, selected with `--ont`.
- **Paired-read input:** two pairs of files supplied through the Hi-C input options `--h1` and `--h2`.
- **Threads:** 128.

In this recorded invocation, all long-read FASTQ files are positional inputs under the `--ont` preset. The separate ultra-long integration option `--ul` is unused. Input filenames are retained from the original command; their names alone do not establish sequencing technology or upstream filtering criteria.

Parameter meanings below were checked against the help output of the same executable. Unspecified settings follow that version and its selected preset.

## Input files

Place the following files in the assembly working directory, or provide links to their actual locations.

| Role in the command | Files |
| --- | --- |
| Hi-C R1 input list, `--h1` | `SRR9026577_1.fq.f.gz`, `SRR9026578_1.fq.f.gz` |
| Hi-C R2 input list, `--h2` | `SRR9026577_2.fq.f.gz`, `SRR9026578_2.fq.f.gz` |
| Positional long-read input 1 | `round2.fq.gz` |
| Positional long-read input 2 | `ul10.pass.fq.gz` |
| Positional long-read inputs 3–11 | `ul1.pass.fq.gz` through `ul9.pass.fq.gz`, in numerical order |
| Positional long-read input 12 | `all.filtered.fq.gz` |

Within the two comma-separated Hi-C lists, keep corresponding R1 and R2 files in the same dataset order. FASTQ inputs are gzip-compressed in the recorded command.

## Recorded command

Make the verified hifiasm executable available as `hifiasm` on `PATH`, or replace `hifiasm` below with its full path. Run from the directory containing the inputs.

```bash
hifiasm \
  -t 128 \
  -l 0 \
  --hom-cov 340 \
  --n-hap 4 \
  --h1 SRR9026577_1.fq.f.gz,SRR9026578_1.fq.f.gz \
  --h2 SRR9026577_2.fq.f.gz,SRR9026578_2.fq.f.gz \
  --ont \
  round2.fq.gz \
  ul10.pass.fq.gz \
  ul1.pass.fq.gz \
  ul2.pass.fq.gz \
  ul3.pass.fq.gz \
  ul4.pass.fq.gz \
  ul5.pass.fq.gz \
  ul6.pass.fq.gz \
  ul7.pass.fq.gz \
  ul8.pass.fq.gz \
  ul9.pass.fq.gz \
  all.filtered.fq.gz \
  --dual-scaf \
  --telo-m CCCTAAA
```

## Explicit parameters

| Parameter | Value | Meaning in hifiasm 0.25.0-r726 |
| --- | --- | --- |
| `-t` | `128` | Number of threads. |
| `-l` | `0` | Disable duplicate purging. Other graph-cleaning operations retain the selected version/preset settings. |
| `--hom-cov` | `340` | Set homozygous read coverage to 340. |
| `--n-hap` | `4` | Set the number of haplotypes to four. |
| `--h1` | Two comma-separated FASTQ files | Supply the Hi-C R1 input list. |
| `--h2` | Two comma-separated FASTQ files | Supply the matching Hi-C R2 input list. |
| `--ont` | Enabled | Select the ONT simplex assembly preset. |
| `--dual-scaf` | Enabled | Output scaffolding. |
| `--telo-m` | `CCCTAAA` | Set the 5′-end telomere-identification motif. |

The haplotype-count setting specifies an assembler parameter; completeness and phasing accuracy require evaluation of the resulting assembly. The telomere option specifies a recognition motif and does not itself establish a telomere-to-telomere result.

## Using the repository script

Run the script inside a PBS job or on a compute node with an appropriate allocation. The command requests 128 threads; request matching CPU resources and choose memory and walltime appropriate to the dataset and cluster policy.

```bash
# Paths are examples to replace for your environment.
# Include this invocation in your PBS job script.
HIFIASM=/path/to/hifiasm-0.25.0/hifiasm \
  bash /path/to/phasing/scripts/01_initial_assembly.sh /path/to/assembly_inputs \
  > /path/to/assembly_inputs/assembly.stdout.log \
  2> /path/to/assembly_inputs/assembly.stderr.log
```

The wrapper changes into the supplied working directory and invokes hifiasm with the recorded parameters. If `HIFIASM` is omitted, it uses `hifiasm` from `PATH`. The two log redirections are a reproducibility convenience in this example; the original command did not specify redirection.

## Outputs and records to retain

The original command does not set `-o`, so the verified executable's default output prefix is **`hifiasm.asm`**. Outputs are written relative to the assembly working directory.

Retain the generated GFA files, standard output/error logs, the exact invocation, and scheduler resource records. Also record the input datasets and their upstream processing. A specific GFA should be identified explicitly when passing the assembly into a downstream graph-processing or phasing step.

This page covers the initial assembly invocation. Subsequent graph pruning, collapsed-node handling, and progeny-based phasing need their own commands and parameter records.

## Software documentation

- [hifiasm source repository](https://github.com/chhylp123/hifiasm)
- [hifiasm documentation](https://hifiasm.readthedocs.io/)

Use the help and source for the recorded version when reproducing this step, since software options can change between releases.
