#!/bin/bash
SAMPLE_FILE=$1
SAMPLE="$( cat "$SAMPLE_FILE" | awk '{print $1}' | cut -d'_' -f1,2)"
CHR="$( cat "$SAMPLE_FILE" | cut -d'_' -f3)"
REF="$( cat "$SAMPLE_FILE" | cut -d'_' -f4,5)"
COMP="$( cat "$SAMPLE_FILE" | cut -d'_' -f6)"
echo $SAMPLE $CHR
echo "REFERENCE" $REF
Rscript --vanilla ~/bin/fruit_DMR_methylKit.R $SAMPLE $CHR $REF $COMP
