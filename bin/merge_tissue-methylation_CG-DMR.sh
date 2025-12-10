TO_LIST="$( cat $1 )"
#echo $TO_LIST
CHR="$( echo $TO_LIST | cut -d'_' -f2 )"
echo "------" $CHR
SAMPLEPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ad_merge_CG-DMR/aa_data"
CHR_SIZE="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/chr.size.bed"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ad_merge_CG-DMR/bb_output"
SAMPLE="$( ls $SAMPLEPATH/*$TO_LIST*.methylation.bed | tr '\n' '\t' )"
KEEP="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/loci"
#echo $SAMPLE
bedtools unionbedg -i $SAMPLE -empty -g $CHR_SIZE -filler NA |\
intersectBed -a - -b $KEEP/$CHR"_CG-DMR-loci_collapsed.bed" > $OUTPATH/$CHR"_CG-DMR.merged.methylation.tmp"

#ls $SAMPLEPATH/*$TO_LIST*.methylation.bed |\
#sed 's./mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ad_merge_CG-DMR/aa_data/..g' |\
#sed 's/_CG_ch01.methylation.bed//g' > $OUTPATH/00_colNames.tsv

cat $OUTPATH/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }'| awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $OUTPATH/$CHR"_CG-DMR.merged.methylation.tmp" |  sed 's/ /\t/g' > $OUTPATH/$CHR"_CG-DMR.merged.methylation.bed"

rm $OUTPATH/$CHR"_CG-DMR.merged.methylation.tmp"
