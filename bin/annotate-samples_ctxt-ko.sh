#!/bin/bash
echo 'Activate gawk environment'
SAMPLE="$( cat $1 | cut -d'_' -f1,2,3,4)" # TS-9_leaf_Biseq_CHH
KO="$( cat $1 | cut -d'_' -f5,6)"
INDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ac_ctxt-conditional-GWAS/ba_phenotype/ca_samples-bed-files"
KODIR="/mnt/disk2/vibanez/04_natural-experimental_DMR/aa_KO"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ac_ctxt-conditional-GWAS/ba_phenotype/cc_samples-with-KOs-features"

echo $SAMPLE "-" $KO
if [[ $ANNO == "general" ]]; then
	intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $KODIR/$KO".ctxt-merged" -wb | gawk '{OFS="\t"}{ if ($18<0) {print}}' > $OUTDIR/$SAMPLE"_"$KO".methylation"
	gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t"  $OUTDIR/$SAMPLE"_"$KO".methylation" > $OUTDIR/$SAMPLE"_"$KO".bed"
else
	intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/${ANNO}.anno |\
	gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_"${ANNO}.bed
fi
