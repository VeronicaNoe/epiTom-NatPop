#!/bin/bash
GROUP="$( cat $1  | cut -d'_' -f1)"
CTXT="$( cat $1 | cut -d'_' -f2 )"
CHR="$( cat $1 | cut -d'_' -f3 )"

echo "------" $CHR
SAMPLEPATH="/mnt/disk3/vibanez/metaplot_tomato_group/data"
CHR_SIZE="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/chr.size.bed"
OUTPATH="/mnt/disk3/vibanez/metaplot_tomato_group/merged_data"
SAMPLE="$( ls $SAMPLEPATH/$GROUP*$CTXT"_"$CHR.bed | tr '\n' '\t' )"

#echo $SAMPLE
bedtools unionbedg -i $SAMPLE -empty -g $CHR_SIZE -filler NA |\
intersectBed -a - -b $SAMPLEPATH/$GROUP"_"$CTXT"_"$CHR".loci_collapsed.bed" > $OUTPATH/$GROUP"_"$CTXT"_"$CHR.methylation.merged.bed

ls $SAMPLEPATH/$GROUP*$CTXT"_"$CHR.bed | sed 's./mnt/disk3/vibanez/metaplot_tomato_group/data..g' |\
cut -d'_' -f2 > $OUTPATH/$GROUP"_"$CTXT"_"$CHR"_colNames.tsv"
