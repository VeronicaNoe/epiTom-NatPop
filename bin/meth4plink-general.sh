#!/bin/bash
SAMPLE="$( echo $1 | sed 's/.merged.methylation.bed//g' )"
echo $SAMPLE
OUTPATH="/mnt/disk2/vibanez/02_methylkit/af_mVCF/aa_data"

cat $SAMPLE".merged.methylation.bed" | sed 's/NA/.|./g' > $OUTPATH/$SAMPLE"_general.meth4plink"
