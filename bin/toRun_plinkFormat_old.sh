#!/bin/bash
SAMPLE_FILE=$1
SAMPLE="$( cat "$SAMPLE_FILE" | cut -d '_' -f1,2 )"
CHR="$( cat "$SAMPLE_FILE" | cut -d '_' -f1)"
ANNO="$( cat "$SAMPLE_FILE" | cut -d '_' -f2)"
DMR="$( cat "$SAMPLE_FILE" | cut -d '_' -f3)"
THRESHOLD="$( cat "$SAMPLE_FILE" | cut -d '_' -f4)"

echo $SAMPLE"."$DMR".meth4plink"
echo $ANNO
echo $CHR
echo $DMR
echo $THRESHOLD
Rscript --vanilla ~/bin/plinkFormat.R $SAMPLE"."$DMR".meth4plink" $ANNO $CHR $DMR $THRESHOLD

