#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
# for DMR-GWAS
Rscript --vanilla ~/bin/check_GIF_DMR-GWAS.R $SAMPLE
