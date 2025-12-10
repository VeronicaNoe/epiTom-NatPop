#!/bin/bash
SAMPLE="$( cat $1  )"
SAMPLE_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.3_methylation-calling"
echo $SAMPLE
cat ~/bin/contexts.tmp | while read line; do
echo $line
zcat ${SAMPLE_DIR}/$SAMPLE.*.CX_report.txt.gz | gawk -v ctxt="$line" '( $6 == ctxt )' > $SAMPLE"_"$line.bed
done
