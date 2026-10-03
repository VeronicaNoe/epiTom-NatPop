#!/bin/bash
SAMPLE="$( cat $1  )"
SAMPLE_DIR="/mnt/ssd123/vibanez/19_ont-mapping/ac_modkit/cx_report"
echo $SAMPLE
cat ~/bin/contexts.tmp | while read line; do
echo $line
zcat ${SAMPLE_DIR}/${SAMPLE}.CX_report.txt.gz | gawk -v ctxt="$line" '( $6 == ctxt )' > $SAMPLE"_"$line.bed
done
