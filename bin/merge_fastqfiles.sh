#!/bin/bash
BARCODE="$( cat $1 | awk '{print $2}' )"
SAMPLE="$( cat $1 | awk '{print $1}' )"
SAMPLE_PATH="/mnt/disk2/vibanez/tomato_methylomes/aa_data"
echo "-------------- Concatenating " $SAMPLE
ls $SAMPLE_PATH/*$BARCODE"-"*"_1.fq.gz"


cat $SAMPLE_PATH/*$BARCODE"-"*"_1.fq.gz" >> $SAMPLE"_1.fq.gz"
cat $SAMPLE_PATH/*$BARCODE"-"*"_2.fq.gz" >> $SAMPLE"_2.fq.gz"
