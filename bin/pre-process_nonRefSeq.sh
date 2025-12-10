#!/bin/bash
SAMPLE="$( echo $1 | sed 's/.filtered.bed//')"
echo $SAMPLE
OUTDIR="/mnt/disk2/vibanez/02_methylkit/nonRefSeq/aa_data"
cat $1 | awk '{OFS="\t"}{print $1, $2, $2, $4/($4+$5)}' | sed 's/,/./g' > $OUTDIR/$SAMPLE

