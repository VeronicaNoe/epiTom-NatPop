#!/bin/bash
OUT_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree"
INDIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.0_get-DMR/ab_predigree"
cat ${INDIR}/*_CG-DMR.bed | sortBed -i - | gawk '{OFS="\t"}{print $1,$2,$3}'|\
mergeBed -i - > ${OUT_DIR}/"CG-DMR-loci_collapsed.bed"

cat ${INDIR}/*_C-DMR.bed | sortBed -i - | gawk '{OFS="\t"}{print $1,$2,$3}'|\
mergeBed -i - > ${OUT_DIR}/"C-DMR-loci_collapsed.bed"
