#!/bin/bash
SAMPLE1="$( cat $1 | cut -d'_' -f1)"
SAMPLE2="$( cat $1 | cut -d'_' -f2)"
CHR="$( cat $1 | cut -d'_' -f3)"
echo  $CHR
Rscript --vanilla ~/bin/pedigree_DMR_methylKit.R $SAMPLE1 $SAMPLE2 $CHR
