#!/bin/bash
SAMPLE="$( cat $1 | cut -d'_' -f1,2,3,4)" # TS-9_leaf_Biseq_CHH_gene
ANNO="$( cat $1 | cut -d'_' -f5)"
INDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ac_ctxt-conditional-GWAS/ba_phenotype/ca_samples-bed-files"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ac_ctxt-conditional-GWAS/ba_phenotype/cb_samples-with-features"

echo ${SAMPLE}
if [[ $ANNO == "general" ]]; then
echo $ANNO
  gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" \
    $INDIR/$SAMPLE.bed.tmp > $OUTDIR/${SAMPLE}_general.bed
else
  echo $ANNO
  intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/${ANNO}.anno | \
    gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" \
    > $OUTDIR/${SAMPLE}_${ANNO}.bed
fi

