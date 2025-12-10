#!/bin/bash
SAMPLE="$( cat $1 | cut -d'_' -f2,3 )"
CHR="$( cat $1 | cut -d'_' -f5 )"
CTXT="$( cat $1 | cut -d'_' -f4 )"
GROUP="$( cat $1  | cut -d'_' -f1 )"
# directory
OUTPATH="/mnt/disk3/vibanez/metaplot_tomato_group/data"
TMPFILES="/mnt/disk3/vibanez/metaplot_tomato_group/tmp"

cat $TMPFILES/$GROUP"_"$SAMPLE"_"$CTXT"_"$CHR".tmp" |\
intersectBed -wb -a - -b $OUTPATH/$GROUP"_"$CTXT"_"$CHR".loci_collapsed.bed" -nonamecheck |\
gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g' > $OUTPATH/$GROUP"_"$SAMPLE"_"$CTXT"_"$CHR".bed"
