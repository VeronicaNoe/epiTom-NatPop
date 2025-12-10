#!/bin/bash
SAMPLE="$( echo $1 | cut -d'_' -f1)"
CHR="$( echo $1 | cut -d'.' -f1 | cut -d'_' -f4 )"
CTXT="$( echo $1 | cut -d'_' -f3 )"
TIS="$( echo $1 | cut -d'_' -f2 )"
METH="$( cat $1 | awk '{print $1}' )"
SD="$( cat $1 | awk '{print $2}' )"
INDIR="/mnt/disk2/vibanez/01_bismark_analysis/04_processingFiles/ac_filter/bb_output"
OUTDIR="/mnt/disk2/vibanez/02_methylkit/FRUITs/descriptive_stat"

echo $SAMPLE
echo $CHR
echo $CTXT
echo $TIS
#create the bed for the DMR
echo -e "$SAMPLE\t$CHR\t$CTXT\t$TIS\t$METH\t$SD" >> $OUTDIR/"leaves-fruit.meanMeth"
