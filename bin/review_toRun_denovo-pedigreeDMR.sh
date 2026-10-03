#!/bin/bash
SAMPLE1="$( cat $1 | cut -d'_' -f1 )"
SAMPLE2="$( cat $1 | cut -d'_' -f2 )"
CTXT="$( cat $1 | cut -d'_' -f3 )"
CTG="$( cat $1 | cut -d'_' -f4 )"
echo  "$SAMPLE1 vs $SAMPLE2"
Rscript --vanilla ~/bin/review_denovo-pedigree_DMR_methylKit.R $SAMPLE1 $SAMPLE2 $CTXT $CTG
