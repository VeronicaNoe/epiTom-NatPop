#!/bin/bash
SAMPLE="$( cat $1 | cut -d'.' -f1 )"
CHR="$( cat $1 | cut -d'_' -f1 | sed 's/ch//g' )"
DMR="$( cat $1 | cut -d'_' -f2 )"
POS="$( cat $1 | cut -d'_' -f3 )"
INDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/sig"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/ba_snp-blocks/cb_result"

echo $SAMPLE
cat $INDIR/$SAMPLE'.mQTL' | tr ':' '\t' | awk -v SAMPLE=$SAMPLE '{OFS="\t"}{print $1, $2, $2+1, $3, $4, $5, SAMPLE}' > $OUTDIR/$SAMPLE".edited.mQTL"
