#!/bin/bash
GROUP="$( cat $1  | cut -d'_' -f1 )"
CTXT="$( cat $1 | cut -d'_' -f3 )"
CHR="$( cat $1 | cut -d'_' -f4 )"
# directory
OUTPATH="/mnt/disk3/vibanez/metaplot_tomato_group/data"
TMPFILES="/mnt/disk3/vibanez/metaplot_tomato_group/tmp"
export TMPDIR=/mnt/disk3/vibanez/tempdir

cat $TMPFILES/$GROUP*$CTXT"_"$CHR".tmp" | sort -k 1.4n -k 2n |\
mergeBed -i - > $OUTPATH/$GROUP"_"$CTXT"_"$CHR".loci_collapsed.bed"
