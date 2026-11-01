#!/usr/bin/env bash
# Supporting commands for genome assessment. Submit compute tasks through PBS.
# Usage: bash 09_genome_quality.sh STAGE ...
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
STAGE="${1:?Specify an evaluation stage}"
shift
case "$STAGE" in
  basic)
    "$PYTHON" "$SCRIPT_DIR/quality/assembly_metrics.py" "$@"
    ;;
  ngs)
    REF="${1:?Assembly FASTA}"; R1="${2:?NGS R1}"; R2="${3:?NGS R2}"
    bwa index "$REF"
    bwa mem -t 40 "$REF" "$R1" "$R2" | samtools sort -@ 10 -o ngs.bam
    samtools index ngs.bam
    samtools flagstat ngs.bam > ngs.flagstat.txt
    samtools coverage ngs.bam > ngs.coverage.tsv
    ;;
  winnowmap)
    REF="${1:?Assembly FASTA}"; READS="${2:?Long reads}"; TYPE="${3:?hifi or ul}"
    case "$TYPE" in hifi) PRESET=map-pb ;; ul) PRESET=map-ont ;; *) exit 1 ;; esac
    meryl count k=15 output repetitive.meryl "$REF"
    meryl print greater-than distinct=0.9998 repetitive.meryl > repetitive_k15.txt
    winnowmap -W repetitive_k15.txt -t 40 -ax "$PRESET" "$REF" "$READS" > "$TYPE.w.sam"
    samtools sort -@ 40 -o "$TYPE.w.bam" "$TYPE.w.sam"
    samtools index "$TYPE.w.bam"
    paftools.js sam2paf "$TYPE.w.sam" > "$TYPE.w.paf"
    ;;
  partition)
    PAF="${1:?Mapper-ranked PAF}"; LIST="${2:?One gzip FASTQ path per line}"
    TYPE="${3:?hifi or ul}"
    "$PYTHON" "$SCRIPT_DIR/quality/assign_haplotypes.py" "$PAF"
    for h in 1 2 3 4; do
      "$PYTHON" "$SCRIPT_DIR/quality/extract_reads_from_list.py" "$LIST" "hap$h.reads.txt" \
        | pigz -c > "hap$h.$TYPE.fq.gz"
    done
    ;;
  remap)
    REF="${1:?Haplotype FASTA}"; READS="${2:?Haplotype-specific reads}"
    TYPE="${3:?hifi or ul}"; PREFIX="${4:?Output prefix}"
    case "$TYPE" in hifi) PRESET=map-hifi ;; ul) PRESET=map-ont ;; *) exit 1 ;; esac
    minimap2 -t 40 -ax "$PRESET" "$REF" "$READS" > "$PREFIX.sam"
    samtools sort -@ 40 -o "$PREFIX.bam" "$PREFIX.sam"
    samtools index "$PREFIX.bam"
    paftools.js sam2paf "$PREFIX.sam" | sort -k6,6V -k8,8n > "$PREFIX.paf"
    ;;
  gci)
    REF="${1:?Haplotype FASTA}"; HAP="${2:?Haplotype number 1-4}"
    HBAM="${3:?HiFi BAM}"; HPAF="${4:?HiFi PAF}"
    UBAM="${5:?UL BAM}"; UPAF="${6:?UL PAF}"
    CHRS="chr1.$HAP,chr2.$HAP,chr3.$HAP,chr4.$HAP,chr5.$HAP,chr6.$HAP,chr7.$HAP,chr8.$HAP"
    GCI -r "$REF" --hifi "$HBAM" "$HPAF" --nano "$UBAM" "$UPAF" \
      --chrs "$CHRS" -t 10 -p -it pdf -mq 0 -d "gci.hap$HAP"
    ;;
  effective)
    REF="${1:?Haplotype FASTA}"; BREAKS="${2:?GCI uncovered-region BED}"
    bedtools maskfasta -fi "$REF" -fo effective.masked.fa -bed "$BREAKS"
    "$PYTHON" "$SCRIPT_DIR/quality/assembly_metrics.py" effective.masked.fa > effective.metrics.json
    ;;
  craq)
    REF="${1:?Haplotype FASTA}"; NGS="${2:?NGS BAM}"; LONG="${3:?Long-read BAM}"
    OUT="${4:?Output directory}"
    craq -ngs "$NGS" -sms "$LONG" -g "$REF" -q 0 -t 10 -D "$OUT"
    ;;
  merqury)
    REF="${1:?Assembly FASTA}"; R1="${2:?NGS R1}"; R2="${3:?NGS R2}"
    : "${MERQURY:?Set MERQURY to the Merqury v1.3 installation directory}"
    export OMP_NUM_THREADS=40
    meryl count k=21 output ngs.k21.meryl "$R1" "$R2"
    bash "$MERQURY/merqury.sh" ngs.k21.meryl "$REF" assembly
    ;;
  busco)
    MODE="${1:?genome or proteins}"; INPUT="${2:?FASTA input}"; OUT="${3:?Output name}"
    : "${BUSCO_LINEAGE:?Set BUSCO_LINEAGE to the local embryophyta_odb10 dataset}"
    busco -m "$MODE" -i "$INPUT" -o "$OUT" -l "$BUSCO_LINEAGE" -c 10 --offline
    ;;
  depth)
    BAM="${1:?Alignment BAM}"; PREFIX="${2:?Output prefix}"
    samtools index "$BAM"
    mosdepth -n --fast-mode --by 100000 -t 40 "$PREFIX" "$BAM"
    ;;
  switch)
    CONTACTS="${1:?MND-format contact records}"; FAI="${2:?Assembly FAI}"; PREFIX="${3:?Output prefix}"
    "$PYTHON" "$SCRIPT_DIR/quality/window_contacts.py" "$CONTACTS" > "$PREFIX.window_pairs.tsv"
    perl "$SCRIPT_DIR/quality/count_window_pairs.pl" "$PREFIX.window_pairs.tsv" > "$PREFIX.cor.tsv"
    "$PYTHON" "$SCRIPT_DIR/quality/sum_haplotype_support.py" "$PREFIX.cor.tsv" > "$PREFIX.stat.tsv"
    awk '$2!=$3 && $4+$5+$6+$7>200' "$PREFIX.stat.tsv" | cut -f 1 \
      | sed 's/_/\t/g' | sort -k1,1 -k2,2n \
      | awk -v OFS="\t" '{print $1,$2*10000,($2+1)*10000}' > "$PREFIX.candidates.bed"
    bedtools merge -i "$PREFIX.candidates.bed" -d 50000 \
      | awk '$3-$2>10000' > "$PREFIX.switch_regions.bed"
    "$PYTHON" "$SCRIPT_DIR/quality/switch_summary.py" "$FAI" "$PREFIX.switch_regions.bed" \
      > "$PREFIX.switch_error_rates.tsv"
    ;;
  *) printf 'Unknown evaluation stage: %s\n' "$STAGE" >&2; exit 1 ;;
esac
