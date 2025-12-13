#!/bin/bash
SAMPLE="$( cat $1 )"
SAMPLE_DIR=$2
OUT_DIR=$3

# intersect windows and samples
echo "------ get C-DMR"
cat ${SAMPLE_DIR}/*${SAMPLE}"_CH"*"DMR.methylation.bed"  | sortBed -i - | mergeBed -i - > ${OUT_DIR}/${SAMPLE}"_C-DMR.tmp"
echo "------ get CG-DMR"
cat ${SAMPLE_DIR}/*${SAMPLE}"_CG_DMR.methylation.bed"  | sortBed -i - | mergeBed -i - |\
intersectBed -a - -b ${OUT_DIR}/${SAMPLE}"_C-DMR.tmp" -v > ${OUT_DIR}/${SAMPLE}"_CG-DMR.tmp"
###
cat ${SAMPLE_DIR}/*${SAMPLE}"_CH"*"DMR.methylation.bed"  |\
intersectBed -a - -b ${OUT_DIR}/${SAMPLE}"_C-DMR.tmp" > ${OUT_DIR}/${SAMPLE}"_C-DMR.bed"
cat ${SAMPLE_DIR}/*${SAMPLE}"_CG_DMR.methylation.bed"  |\
intersectBed -a - -b ${OUT_DIR}/${SAMPLE}"_CG-DMR.tmp" > ${OUT_DIR}/${SAMPLE}"_CG-DMR.bed"

rm ${OUT_DIR}/${SAMPLE}*tmp
