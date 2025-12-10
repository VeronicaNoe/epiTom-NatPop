#!/bin/bash
echo 'Activate gawk environment'
KODIR="/mnt/disk2/vibanez/04_natural-experimental_DMR/aa_KO"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTDIR="/mnt/disk2/vibanez/04_natural-experimental_DMR/annotate_kos"

KO="$( cat $1 | cut -d'_' -f1)"
ANNO="$( cat $1 | cut -d'_' -f2)"

echo $KO "-" $ANNO
intersectBed -a $KODIR/$KO.ctxt-merged.hypomethylated -b $ANNODIR/$ANNO.merged.anno  > $OUTDIR/$KO"_"$ANNO.bed

#for Old
#KO="$( cat $1 | cut -d'_' -f1,2)"
#ANNO="$( cat $1 | cut -d'_' -f3)"
#FEATKO="$( cat $1 | cut -d'_' -f4 | cut -d'.' -f1 )"
#CTXT="$( cat $1 | cut -d'_' -f2 | cut -d'-' -f1 )"

#echo $KO "-" $ANNO "-" $FEATKO ""
#gawk '{OFS="\t"}{if ($14<0) {print}}' $KODIR/$KO.methylation.bed |\
#intersectBed -a - -b $ANNODIR/$ANNO.anno | intersectBed -a - -b $KODIR/$KO.methylation.bed |\
#gawk '{sum5 += $5; sum6 += $6; sum8 += $8; sum9 += $9} END {print sum5, sum6, sum6/sum5, sum8, sum9, sum9/sum8, $11}' OFS="\t" > $OUTDIR/$KO"_"$ANNO"_"$FEATKO.bed
