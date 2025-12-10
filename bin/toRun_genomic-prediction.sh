#!/bin/bash
MM="$( cat "$1" | cut -d'_' -f1)"
TRAIT="$( cat "$1" | cut -d '_' -f2)"
TISSUE="$( cat "$1" | cut -d '_' -f3)"
DATASET="$( cat "$1" | cut -d'_' -f4)"

### This is old, check and remove it #Rscript --vanilla /mnt/disk2/vibanez/10_data-analysis/Genomic_prediction/01.2_Tomato_models.R $MM $TRAIT $TISSUE
# WORKINKG -> Rscript --vanilla ~/bin/fig5_02.3_predictive-models.R $MM $TRAIT $TISSUE $DATASET
Rscript --vanilla ~/bin/run4ptnt_rrreml_loco.R $MM $TRAIT $TISSUE $DATASET # only for ptnt
