#!/bin/bash
SAMPLE="$( cat $1 )"
SAMPLE_DIR=$2
#"/mnt/disk2/vibanez/04_methylome-comparison"
OUT_DIR=$3
COV_DIR="cov3"

#"/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.0_get-DMR/ab_output"
# intersect windows and samples
echo "------ get C-DMR"
cat ${SAMPLE_DIR}/${COV_DIR}/*${SAMPLE}"_CH"*"DMR.methylation.bed"  | sortBed -i - | mergeBed -i - > ${OUT_DIR}/${SAMPLE}"_C-DMR.tmp"
echo "------ get CG-DMR"
cat ${SAMPLE_DIR}/${COV_DIR}/*${SAMPLE}"_CG_DMR.methylation.bed"  | sortBed -i - | mergeBed -i - |\
intersectBed -a - -b ${OUT_DIR}/${SAMPLE}"_C-DMR.tmp" -v > ${OUT_DIR}/${SAMPLE}"_CG-DMR.tmp"
###
cat ${SAMPLE_DIR}/${COV_DIR}/*${SAMPLE}"_CH"*"DMR.methylation.bed"  |\
intersectBed -a - -b ${OUT_DIR}/${SAMPLE}"_C-DMR.tmp" > ${OUT_DIR}/${SAMPLE}"_C-DMR.bed"
cat ${SAMPLE_DIR}/${COV_DIR}/*${SAMPLE}"_CG_DMR.methylation.bed"  |\
intersectBed -a - -b ${OUT_DIR}/${SAMPLE}"_CG-DMR.tmp" > ${OUT_DIR}/${SAMPLE}"_CG-DMR.bed"

rm ${OUT_DIR}/${SAMPLE}*tmp
