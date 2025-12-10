#!/bin/bash
SAMPLE="$( cat $1 )"
DMR="$( cat $1 | cut -d'_' -f2)"
CHR="$( cat $1 | cut -d'_' -f1)"

echo $SAMPLE
Rscript --vanilla ~/bin/get_meth-expression_correlations.R $CHR $DMR $SAMPLE
