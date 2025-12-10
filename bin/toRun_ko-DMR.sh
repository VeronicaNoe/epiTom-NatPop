#!/bin/bash
SAMPLE_FILE=$1
SAMPLE="$( cat "$SAMPLE_FILE" | cut -d'_' -f1,2,3)"
CHR="$( cat "$SAMPLE_FILE" | cut -d'_' -f4)"
echo $SAMPLE $CHR
Rscript --vanilla ~/bin/ko_DMR_methylkit.R $SAMPLE $CHR
