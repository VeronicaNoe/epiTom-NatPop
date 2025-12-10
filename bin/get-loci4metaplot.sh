#!/bin/bash
SAMPLE="$( cat $1 | cut -d'_' -f2,3 )"
CHRCTXT="$( cat $1 | cut -d'_' -f4,5)"
GROUP="$( cat $1  | cut -d'_' -f1 )"
echo $CHRCTXT
# directory
SAMPLEPATH="/mnt/disk2/vibanez/01_bismark_analysis/04_processingFiles/ac_filter/bb_output"
TMPFILES="/mnt/disk3/vibanez/metaplot_tomato_group/tmp"

# creating a small data set with all the windows
cat $SAMPLEPATH/$SAMPLE"_Biseq_"$CHRCTXT".filtered.bed" | gawk '{OFS="\t"}{print $1,$2,$2+1,$4,($4+$5)}' |\
sed 's/,/./g' > $TMPFILES/$GROUP"_"$SAMPLE"_"$CHRCTXT.tmp
