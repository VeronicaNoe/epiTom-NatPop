#!/bin/bash
SAMPLE="$( cat "$1" )"
SAMPLE_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.4_filtering/ab_chr-split/bb_data"

cat ${SAMPLE_DIR}/$SAMPLE.bed | grep -vE "SL2.50ch01|SL2.50ch02|SL2.50ch03|SL2.50ch04|SL2.50ch05|SL2.50ch06|SL2.50ch07|SL2.50ch08|SL2.50ch09|SL2.50ch10|SL2.50ch11|SL2.50ch12" > $SAMPLE"_chrNonRef.bed"
cat ~/bin/chromosomes.tmp | while read line; do
	echo $line
	cat ${SAMPLE_DIR}/${SAMPLE}.prefiltered.bed | grep -e $line > $SAMPLE"_"$line".bed"
done
