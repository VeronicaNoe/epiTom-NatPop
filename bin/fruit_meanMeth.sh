#!/bin/bash
SAMPLE="$( echo $1 | cut -d'_' -f1)"
CHR="$( echo $1 | cut -d'.' -f1 | cut -d'_' -f5 )"
CTXT="$( echo $1 | cut -d'_' -f4 )"
TIS="$( echo $1 | cut -d'_' -f2 )"

INDIR="/mnt/disk2/vibanez/01_bismark_analysis/04_processingFiles/ac_filter/bb_output"
OUTDIR="/mnt/disk2/vibanez/02_methylkit/FRUITs/descriptive_stat/aa_data"

echo $SAMPLE
echo $CHR
echo $CTXT
echo $TIS
#create the bed for the DMR
cat $INDIR/$1 | awk '{print $4/($4+$5)}' | datamash mean 1 sstdev 1 | tr ',' '.' > $OUTDIR/$SAMPLE"_"$CHR"_"$CTXT"_"$TIS".meanMeth"
