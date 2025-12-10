#!/bin/bash
GENE="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/allGenes.bed"
KO_DIR="/mnt/disk2/vibanez/09_KO-processing/09.0_biseq/ae_DMR/ba_DMR-calling"
SAMPLE="$( cat "$1" )"
OUT_DIR=$2

sortBed -i $GENE | closestBed -a ${KO_DIR}/$SAMPLE.ctxt-merged -b - -D a |\
 awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $18, $19}' > ${OUT_DIR}/${SAMPLE}.closest-gene
