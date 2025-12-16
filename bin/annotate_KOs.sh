#!/bin/bash
KODIR="09.3_DMR-comparison"
ANNODIR="05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUTDIR="09.3_DMR-comparison/aa_annotate-kos"

KO="$( cat $1 | cut -d'_' -f1)"
ANNO="$( cat $1 | cut -d'_' -f2)"

echo $KO "-" $ANNO
intersectBed -a $KODIR/$KO.ctxt-merged.hypomethylated -b $ANNODIR/$ANNO.merged.anno  > $OUTDIR/$KO"_"$ANNO.bed

