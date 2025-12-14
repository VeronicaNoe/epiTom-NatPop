#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
SAMPLE_PATH="08_fruit-processing/08.0_get-DMR"
OUT="08_fruit-processing/08.1_DMR-classification"
# intersect windows and samples
echo "------ subtract CHG and CHH"
cat $SAMPLE_PATH/$SAMPLE"_CH"*"_DMR.methylation.bed" | sortBed -i - | mergeBed -i - |\
bedtools subtract -a $SAMPLE_PATH/$SAMPLE"_CG_DMR.methylation.bed" -b - > $OUT/$SAMPLE.CG-DMR.bed
