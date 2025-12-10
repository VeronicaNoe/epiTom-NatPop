#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
# --- Define samples to process
SAMPLE_PATH="/mnt/disk2/vibanez/08_fruit-processing/08.0_get-DMR"
TMP="/mnt/disk2/vibanez/08_fruit-processing/08.0_get-DMR/tmp"
OUT="/mnt/disk2/vibanez/08_fruit-processing/08.1_DMR-classification/ab_output"
# intersect windows and samples
echo $SAMPLE
echo "------ subtract CG"
	bedtools subtract -a $SAMPLE_PATH/$SAMPLE"_CG_DMR.methylation.bed" -b $OUT/$SAMPLE".CG-DMR.bed" >> $TMP/$SAMPLE.C-DMR.tmp
echo "------ subtract CHG"
        bedtools subtract -a $SAMPLE_PATH/$SAMPLE"_CHG_DMR.methylation.bed" -b $OUT/$SAMPLE".CG-DMR.bed" >> $TMP/$SAMPLE.C-DMR.tmp
echo "------ subtract CHH"
        bedtools subtract -a $SAMPLE_PATH/$SAMPLE"_CHH_DMR.methylation.bed" -b $OUT/$SAMPLE".CG-DMR.bed" >> $TMP/$SAMPLE.C-DMR.tmp
echo "------ sort"
	bedtools sort -i $TMP/$SAMPLE.C-DMR.tmp > $TMP/$SAMPLE.C-DMR.bed
echo "------ sum with R"
	Rscript --vanilla ~/bin/sum_fruit_windows.R $SAMPLE.C-DMR.bed
