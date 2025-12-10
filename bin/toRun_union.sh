#!/bin/bash
SAMPLE_FILE=$1
SAMPLE="$( cat "$SAMPLE_FILE" )"
echo $SAMPLE
Rscript --vanilla ~/bin/unite_with_methylkit.R $SAMPLE

