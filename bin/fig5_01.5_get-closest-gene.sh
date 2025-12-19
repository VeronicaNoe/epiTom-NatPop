#!/bin/bash
QTL="10_data-analysis/Fig5/results/01.1_QTLs_with_coordinates.tsv"
ANNO_DIR="05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUT_DIR="10_data-analysis/Fig5/results"

# Fix chromosome naming with leading zeros
awk 'NR > 1 {
  OFS = "\t";
  chr = $7;
  if (chr < 10) {
    chr = "0" chr;
  }
  print chr, $8, $9, $1, $2, $3, $4, $5, $6;
}' "${QTL}" |\
sortBed -i - > ${OUT_DIR}/tmp.bed

closestBed -D b -a ${OUT_DIR}/tmp.bed  -b ${ANNO_DIR}/allGeneNames-Function.bed -t all > ${OUT_DIR}/04.0_QTL_closest-gene.bed

intersectBed -a ${OUT_DIR}/tmp.bed  -b ${ANNO_DIR}/TE-wo-gene.bed  -wo |\
closestBed -D b -a -  -b ${ANNO_DIR}/allGeneNames-Function.bed -t all > ${OUT_DIR}/04.1_QTL-overTE_closest-gene.bed

awk '{print $15}' ${OUT_DIR}/04.0_QTL_closest-gene.bed | sed 's/Parent=gene://g' > ${OUT_DIR}/04.2_candidate-gene_list

bedtools window -a ${OUT_DIR}/tmp.bed  -b ${ANNO_DIR}/allGeneNames-Function.bed -w 1000000 > ${OUT_DIR}/04.3_nearby-genes.tsv

while read line; do
	grep -B5 -A 5 "$line" ${ANNO_DIR}/allGeneNames-Function.bed |\
	awk '{OFS="\t"}{print $6, $7}' | sed 's/Parent=gene://g' |\
	awk -v var="$line" '{print var "\t" $0}'
done < ${OUT_DIR}/04.2_candidate-gene_list > ${OUT_DIR}/04.4_nearby-genes2candidate-gene.tsv
