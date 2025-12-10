#!/bin/bash
GENE="$( cat "$1" )"
echo $GENE
Rscript --vanilla ~/bin/fig4_00.2_get_EpiallelicStatus.R $GENE
