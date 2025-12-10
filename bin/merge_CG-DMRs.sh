CHR="$( cat $1 )"
echo "------" $CHR
SAMPLEPATH="/mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/aa_data"
CHR_SIZE="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/chr.size.bed"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/bb_output"
KEEP="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ae_loci_CG-DMRs/tmp"
KEEPC="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/af_loci_C-DMRs/tmp"
CG="$( ls $SAMPLEPATH/$CHR*CG-DMR.bed | tr '\n' '\t' )"
C="$( ls $SAMPLEPATH/$CHR*C-DMR.bed | tr '\n' '\t' )"

echo "======== Procesing CG-DMRs ========"
bedtools unionbedg -i $CG -empty -g $CHR_SIZE -filler NA |\
	intersectBed -a - -b $KEEP/$CHR"_CG-DMR-loci_collapsed.bed" > $OUTPATH/$CHR"_CG-DMR.merged.tmp"

#ls $SAMPLEPATH/$CHR*.bed | sed 's./mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/aa_data/..g' |\
#	sed 's/.CG-DMR.bed//g' > $OUTPATH/"CG_colNames.tsv"

cat $OUTPATH/"CG_colNames.tsv" | awk 'BEGIN { ORS = "\t" } { print }'|\
	awk '{OFS="\t"}{print "chr","start","end", $0}'|\
	cat - $OUTPATH/$CHR"_CG-DMR.merged.tmp" > $OUTPATH/$CHR"_CG-DMR.merged"

rm $OUTPATH/$CHR"_CG-DMR.merged.tmp"

echo "======== Procesing C-DMRs ========"
bedtools unionbedg -i $C -empty -g $CHR_SIZE -filler NA |\
        intersectBed -a - -b $KEEPC/$CHR"_C-DMR-loci_collapsed.bed" > $OUTPATH/$CHR"_C-DMR.merged.tmp"

#ls $SAMPLEPATH/$CHR*C-DMR.bed | sed 's./mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR/aa_data/..g' |\ 
#        sed 's/.C-DMR.bed//g' > $OUTPATH/"C_colNames.tsv"

cat $OUTPATH/"C_colNames.tsv" | awk 'BEGIN { ORS = "\t" } { print }'|\
        awk '{OFS="\t"}{print "chr","start","end", $0}'|\
        cat - $OUTPATH/$CHR"_C-DMR.merged.tmp" > $OUTPATH/$CHR"_C-DMR.merged"

rm $OUTPATH/$CHR"_C-DMR.merged.tmp"
