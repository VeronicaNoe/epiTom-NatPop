#!/bin/bash
MM="$( cat "$1" | cut -d'_' -f1)"
TRAIT="$( cat "$1" | cut -d '_' -f2)"
TISSUE="$( cat "$1" | cut -d '_' -f3)"
DATASET="$( cat "$1" | cut -d'_' -f4)"

Rscript --vanilla ~/bin/fig5_02.3_predictive-models.R $MM $TRAIT $TISSUE $DATASET
