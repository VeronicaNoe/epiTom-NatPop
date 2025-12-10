#!/bin/bash
#!/bin/bash
QTL="/mnt/disk2/vibanez/10_data-analysis/Fig5/results/Tm159_mQTL"
ANNO_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUT_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/results"

# Fix chromosome naming with leading zeros
awk '{
  OFS = "\t";
  chr = $7;
  if (chr < 10) {
    chr = "0" chr;
  }
  print chr, $8, $9, $1, $2, $3, $4, $5, $6;
}' "${QTL}" |\
sortBed -i - > ${OUT_DIR}/tmp_candidate.bed
closestBed -D b -a ${OUT_DIR}/tmp_candidate.bed  -b ${ANNO_DIR}/allGeneNames-Function.bed -t all > ${OUT_DIR}/10.0_candidate_gene_closest-gene.bed
intersectBed -a ${OUT_DIR}/tmp_candidate.bed  -b ${ANNO_DIR}/TE-wo-gene.bed  -wo |\
closestBed -D b -a -  -b ${ANNO_DIR}/allGeneNames-Function.bed -t all > ${OUT_DIR}/10.0_candidate_QTL-overTE_closest-gene.bed
