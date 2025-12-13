#!/bin/bash
SAMPLE="$( cat $1 )"
SAMPLE_DIR=$2
OUT_DIR=$3

# intersect windows and samples
#echo $SAMPLE
cat ${SAMPLE_DIR}/*"-DMR_merged.bed" | sortBed -i - | mergeBed -i - > ${OUT_DIR}/$SAMPLE.C-DMR.bed
