#!/bin/bash
SAMPLE_FILE=$1
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ae_merge_CG-DMRs/aa_data"

cat $SAMPLE_FILE | awk '{ print $1, $2, $4, $9, $10, $11}' | sed 's/ /\t/g' > $OUTPATH/$SAMPLE_FILE
