#!/bin/bash
SAMPLE="$( cat $1 )"
SAMPLE_DIR=$2
OUT_DIR=$3
#TMP="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/tmp-files"
#ls $TMP

# intersect windows and samples
#echo $SAMPLE
cat ${SAMPLE_DIR}/*"-DMR_merged.bed" | sortBed -i - | mergeBed -i - > ${OUT_DIR}/$SAMPLE.C-DMR.bed
