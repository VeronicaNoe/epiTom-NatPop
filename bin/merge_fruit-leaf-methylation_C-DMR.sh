#!/bin/bash
CHR="$( cat $1 | cut -d'_' -f1 )"
echo "------" $CHR
SAMPLEPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ah_merge_C-DMR/bd_output"
CHR_SIZE="/mnt/disk2/vibanez/02_methylkit/FRUITs/ah_merge_C-DMR/chr.size.bed"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ca_dendrogram/C-merged-data/bb_output"
SAMPLE="$( ls $SAMPLEPATH/$CHR*.weighted.methylation.bed | tr '\n' '\t' )"
KEEP="/mnt/disk2/vibanez/02_methylkit/FRUITs/ae_loci_C-DMRs/tmp"
#echo $SAMPLE

bedtools unionbedg -i $SAMPLE -empty -g $CHR_SIZE -filler NA |\
intersectBed -a - -b $KEEP/"fl_"$CHR"_C-DMR-loci_collapsed.bed" > $OUTPATH/$CHR"_C-DMR.merged.methylation.tmp"

ls $SAMPLEPATH/$CHR*.weighted.methylation.bed |\
sed 's./mnt/disk2/vibanez/02_methylkit/FRUITs/ah_merge_C-DMR/bd_output/..g' |\
sed 's/.weighted.methylation.bed//g' | cut -d '_' -f2,3  > $OUTPATH/C_colNames.tsv

cat $OUTPATH/C_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }'|\
awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $OUTPATH/$CHR"_C-DMR.merged.methylation.tmp" > $OUTPATH/$CHR"_C-DMR.merged.methylation.bed"

rm $OUTPATH/$CHR"_C-DMR.merged.methylation.tmp"

