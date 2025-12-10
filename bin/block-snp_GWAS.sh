#!/bin/bash
SAMPLE="$( cat $1 )"
CHR="$( cat $1 | cut -d'_' -f1 | sed 's/ch//g' )"
DMR="$( cat $1 | cut -d'_' -f2 )"
POS="$( cat $1 | cut -d'.' -f1 | cut -d'_' -f3 )"
INDIR="/mnt/disk2/vibanez/GWAS/analysis_Result/sigSNPs"
CISDIR="/mnt/disk2/vibanez/GWAS/analysis_Result/bd_cis-trans"
OUTDIR="/mnt/disk2/vibanez/GWAS/analysis_Result/bc_snp-blocks/bc_result"
DISDIR="/mnt/disk2/vibanez/GWAS/analysis_Result/bh_distance/bc_result"

echo $SAMPLE
#create the bed for the DMR
	touch $DMR'_'$CHR'_'$POS".bed"
	echo $CHR":"$POS":"$POS | tr ':' '\t' | awk '{OFS="\t"}{print $1, $2, $3+99}'> $DMR'_'$CHR'_'$POS".bed"
# get the snps in cis and total
	cat $INDIR/$SAMPLE'.mQTL' | tr ':' '\t' | awk '{OFS="\t"}{print $1, $2, $2+1, $3, $4, $5}' |\
	bedtools window -a $DMR'_'$CHR'_'$POS".bed" -b - -w 100000 > $CISDIR/$DMR'_'$CHR'_'$POS".cis"
        cat $CISDIR/$DMR'_'$CHR'_'$POS".cis" | awk '{OFS="\t"}{print $6-$3}' >  $DISDIR/$DMR'_'$CHR'_'$POS'.distance'
        cat $INDIR/$SAMPLE'.mQTL' | wc -l > $CISDIR/$DMR'_'$CHR'_'$POS'.total'
#make the blocks
	cat $INDIR/$SAMPLE'.mQTL' | tr ':' '\t' | awk '{OFS="\t"}{print $1, $2, $2+1, $3, $4, $5}' | mergeBed -i - -d 1000000 -c 4,5,6 -o mean,mean,min > $OUTDIR/$SAMPLE'.mQTL'

rm $DMR'_'$CHR'_'$POS".bed"
