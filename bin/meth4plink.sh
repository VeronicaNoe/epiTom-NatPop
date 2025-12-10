#!/bin/bash
SAMPLE="$( echo $1 | sed 's/.methylation//g' )"
echo $SAMPLE
OUTPATH="/mnt/disk2/vibanez/02_methylkit/af_mVCF/aa_data"

cat $SAMPLE".methylation" | sed 's/NA/.|./g' > $SAMPLE."meth4plink"  #$OUTPATH/$SAMPLE".meth4plink"
