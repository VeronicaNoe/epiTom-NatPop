#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
# for DMR-GWAS
Rscript --vanilla ~/bin/avg-meth-4-groups-metaplot.R $SAMPLE
