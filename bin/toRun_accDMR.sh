#!/bin/bash
SAMPLE="$( cat $1 |  cut -d'_' -f1,2 )"
CHR="$( cat $1 | cut -d'_' -f3)"
echo $SAMPLE $CHR
Rscript --vanilla bin/acc_DMR_methylKit.R $SAMPLE $CHR
