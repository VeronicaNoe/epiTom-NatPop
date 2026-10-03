#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
# for DMR-GWAS
Rscript --vanilla ~/bin/check_GIF_DMR-GWAS_TIP.R $SAMPLE
mv ca_targets/${SAMPLE}.target cd_targets-done/
