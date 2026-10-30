#!/usr/bin/env bash
# Submit compute-intensive stages through PBS. FASTK v1.2.
# Usage: bash 04_progeny_kmers.sh {parent|sample|merge|markers|contacts|anchor} ...
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
STAGE="${1:?Specify a stage}"
shift

case "$STAGE" in
  parent)
    # Combined parental assembly FASTA; identifiers must be unique across parents.
    PARENTS="${1:?Specify the combined parental FASTA}"
    mkdir -p tmp
    FastK -v -k51 -T10 -M20 -t1 -P./tmp/ "$PARENTS" -NG51
    Logex -T10 'G51uk=A[1]' G51
    Tabex G51uk | awk '$4==1 {print ">kmer"$1"\n"toupper($2)}' > unique_51mers.fa
    "$PYTHON" "$SCRIPT_DIR/progeny/localize_markers.py" "$PARENTS" unique_51mers.fa \
      > parental_markers.fa
    ;;
  sample)
    BASE="${1:?Specify the progeny basename}"
    CATALOG="${2:?Specify the G51uk table prefix}"
    mkdir -p tmp
    FastK -v -k51 -t1 -T10 -M20 -P./tmp/ "${BASE}_R1.f.fq.gz" -N"${BASE}_1"
    FastK -v -k51 -t1 -T10 -M20 -P./tmp/ "${BASE}_R2.f.fq.gz" -N"${BASE}_2"
    # Sum R1/R2 counts, then retain counts >=3 (frequency >2).
    Logex -h1000 -T8 "${BASE}=(A|+B)[3-]" "${BASE}_1" "${BASE}_2"
    Logex -T10 "${BASE}_uG=A&.B" "$BASE" "$CATALOG"
    Tabex "${BASE}_uG" | pigz -p 10 -c > "${BASE}_uG.list.gz"
    ;;
  merge)
    LIST="${1:?Specify the ordered progeny list}"
    # Partition by all 64 trinucleotide prefixes, retaining one common header.
    HEADER_WRITTEN=0
    for a in a c g t; do
      for b in a c g t; do
        for c in a c g t; do
          PART="${a}${b}${c}.merged.txt"
          "$PYTHON" "$SCRIPT_DIR/progeny/merge_vectors.py" "$LIST" "${a}${b}${c}" > "$PART"
          if [[ "$HEADER_WRITTEN" == 0 ]]; then
            head -n 1 "$PART"
            HEADER_WRITTEN=1
          fi
          tail -n +2 "$PART"
        done
      done
    done > progeny.merged.tsv
    ;;
  markers)
    ASSEMBLY="${1:?Specify the combined parental FASTA}"
    MARKERS="${2:?Specify parental_markers.fa}"
    MATRIX="${3:?Specify progeny.merged.tsv}"
    samtools faidx "$ASSEMBLY"
    bedtools makewindows -w 10000 -g "$ASSEMBLY.fai" > parental.10000.bed
    "$PYTHON" "$SCRIPT_DIR/progeny/stat_kmer.py" "$MARKERS" "$MATRIX" > marker_counts.tsv
    "$PYTHON" "$SCRIPT_DIR/progeny/select_window_markers.py" \
      parental.10000.bed marker_counts.tsv 10000 > selected_markers.txt
    "$PYTHON" "$SCRIPT_DIR/progeny/filter_matrix.py" selected_markers.txt "$MATRIX" > selected_vectors.tsv
    "$PYTHON" "$SCRIPT_DIR/progeny/change_kmerid.py" "$MARKERS" selected_vectors.tsv \
      > positioned_vectors.tsv
    ;;
  contacts)
    QUERY="${1:?Specify query vectors with a header}"
    REFERENCE="${2:?Specify reference vectors with a header}"
    "$PYTHON" "$SCRIPT_DIR/progeny/make_pseudocontacts.py" "$QUERY" "$REFERENCE" > progeny.3ddna.txt
    ;;
  anchor)
    GROUP="${1:?Specify the chromosome-group basename}"
    NODES="${2:?Specify the chromosome-group node list}"
    CONTACTS="${3:?Specify progeny.3ddna.txt}"
    UNITIGS="${4:?Specify the unitig FASTA}"
    : "${THREEDDNA:?Set THREEDDNA to the 3D-DNA 201008 installation directory}"
    CONTACTS_ABS="$(realpath -- "$CONTACTS")"
    "$PYTHON" "$SCRIPT_DIR/progeny/filter_group_contacts.py" "$NODES" "$CONTACTS_ABS" > "$GROUP.3ddna.txt"
    seqkit grep -f "$NODES" "$UNITIGS" > "$GROUP.fa"
    FASTA_ABS="$(realpath -- "$GROUP.fa")"
    GROUP_CONTACTS_ABS="$(realpath -- "$GROUP.3ddna.txt")"
    mkdir -p "$GROUP"
    (cd -- "$GROUP"; bash "$THREEDDNA/run-asm-pipeline.sh" -r 1 -i 100 "$FASTA_ABS" "$GROUP_CONTACTS_ABS")
    ;;
  *) printf 'Unknown stage: %s\n' "$STAGE" >&2; exit 1 ;;
esac
