#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
SAMPLE_PATH="/mnt/disk2/vibanez/08_fruit-processing/08.0_get-DMR"
OUT="/mnt/disk2/vibanez/08_fruit-processing/08.1_DMR-classification/ab_output"
# intersect windows and samples
echo "------ subtract CHG and CHH"
cat $SAMPLE_PATH/$SAMPLE"_CH"*"_DMR.methylation.bed" | sortBed -i - | mergeBed -i - |\
bedtools subtract -a $SAMPLE_PATH/$SAMPLE"_CG_DMR.methylation.bed" -b - > $OUT/$SAMPLE.CG-DMR.bed
