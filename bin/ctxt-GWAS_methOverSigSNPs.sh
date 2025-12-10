#!/bin/bash
echo 'Activate gawk environment'
SAMPLE="$( cat $1 )" # TS-9_leaf_Biseq_CHH
DIR="/mnt/disk2/vibanez/descriptive-stat/ac_ctxt-GWAS/annotation"
INDIR="/mnt/disk2/vibanez/descriptive-stat/aa_tmp-bed"
echo $SAMPLE
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $DIR/00_coordinates_genesAssociated.bed |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5), $9, $10, $11}' OFS="\t" > $DIR/$SAMPLE"_methOverAssociatedGene.bed"

intersectBed -b $INDIR/$SAMPLE.bed.tmp -a $DIR/coordinates_genesAssociated.bed -wo |\
gawk '{OFS="\t"} {print $1,$2,$3,$4,$5,$9,$10}' |\
sortBed -i - |\
mergeBed -i - -c 4,5,6,7 -o distinct,distinct,sum,sum > $DIR/methOverAssociated-gene/$SAMPLE"_methOverAssociatedGene.bed"
