#!/bin/bash
GROUP="$( cat $1 | cut -d'_' -f1 )"
DMR="$( cat $1 | cut -d'_' -f2 )"
FEATURE="$( cat $1 | cut -d'_' -f3)"
INDIR="/mnt/data6/vibanez/mVCF/aa_data"
OUTDIR="/mnt/data6/vibanez/mVCF/ab_split-mvcf-group/aa_data"
echo $GROUP
echo $SAMPLE

vcftools --gzvcf $INDIR/$FEATURE.$DMR".mvcf.gz" \
	--keep $INDIR/$GROUP".group" \
	--recode \
	--recode-INFO-all \
	--out $OUTDIR/$GROUP"_"$DMR"_"$FEATURE"_allAcc"

vcftools --gzvcf $INDIR/$FEATURE.$DMR".mvcf.gz" \
        --keep $INDIR/$GROUP".20acc.group" \
        --recode \
        --recode-INFO-all \
        --out $OUTDIR/$GROUP"_"$DMR"_"$FEATURE"_20acc"
