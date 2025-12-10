ANNO="$( cat $1 | cut -d'_' -f1 )"
echo $ANNO
CHR="$( cat $1 | cut -d'_' -f2 )"
echo $CHR
SAMPLEPATH="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/bc_CG-DMR_intersect"
#echo $SAMPLEPATH
#echo $OUTPATH

cat $SAMPLEPATH/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/$ANNO".anno" > $OUTPATH/$CHR"_"$ANNO.CG-DMR.methylation
