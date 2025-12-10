#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
# for DMR-GWAS
Rscript --vanilla ~/bin/combine-manhattan_DMR-GWAS.R $SAMPLE
