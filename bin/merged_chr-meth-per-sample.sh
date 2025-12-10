#!/bin/bash
SAMPLE="$( cat $1 )" # TS-9_leaf_Biseq_CHH_ch*.filtered.bed
INDIR="/mnt/disk2/vibanez/03_biseq-processing/03.4_filtering/ac_filter"
ANNODIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ac_ctxt-conditional-GWAS/ba_phenotype/ca_samples-bed-files"

echo $SAMPLE
cat $INDIR/$SAMPLE"_ch"*".filtered.bed" | gawk '{OFS="\t"} {print $1, $2, $2+1, $4, $5}' |\
sortBed -i - | mergeBed -i - -c 4,5 -o sum,sum > $OUTDIR/$SAMPLE.bed.tmp
