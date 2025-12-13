#!/bin/bash
SAMPLE="$( cat "$1" )"
SAMPLE_DIR="03_biseq-processing/03.4_filtering/ab_chr-split/bb_data"

cat ~/bin/chromosomes.tmp | while read line; do
	echo $line
	cat ${SAMPLE_DIR}/${SAMPLE}.prefiltered.bed | grep -e $line > $SAMPLE"_"$line".bed"
done
