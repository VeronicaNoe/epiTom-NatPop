#!/bin/bash
GROUP="$( cat $1 | cut -d'_' -f1 )"
DMR="$( cat $1 | cut -d'_' -f2 )"
FEATURE="$( cat $1 | cut -d'_' -f3)"
INDIR="06_get-meth-vcf/aa_output"
OUTDIR="06_get-meth-vcf/ad_split-mvcf-group"
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
