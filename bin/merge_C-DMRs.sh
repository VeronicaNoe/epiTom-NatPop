CHR="$( cat $1 )"
echo "------" $CHR
SAMPLEPATH="/mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/aa_data"
CHR_SIZE="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/chr.size.bed"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/bb_output"
KEEPC="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/af_loci_C-DMRs/tmp"
C="$( ls $SAMPLEPATH/$CHR*C-DMR.bed | tr '\n' '\t' )"

echo "======== Procesing C-DMRs ========"
bedtools unionbedg -i $C -empty -g $CHR_SIZE -filler NA |head |\
        intersectBed -a - -b $KEEPC/$CHR"_C-DMR-loci_collapsed.bed" > $OUTPATH/$CHR"_C-DMR.merged.tmp"

#ls $SAMPLEPATH/$CHR*C-DMR.bed | sed 's./mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/aa_data/..g' |\
#        sed 's/.C-DMR.bed//g' > $OUTPATH/"C_colNames.tsv"

cat $OUTPATH/"C_colNames.tsv" | awk 'BEGIN { ORS = "\t" } { print }'|\
        awk '{OFS="\t"}{print "chr","start","end", $0}'|\
        cat - $OUTPATH/$CHR"_C-DMR.merged.tmp" > $OUTPATH/$CHR"_C-DMR.merged"

rm $OUTPATH/$CHR"_C-DMR.merged.tmp"

