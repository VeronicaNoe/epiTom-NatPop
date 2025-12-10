#!/bin/bash
echo 'Activate gawk environment'
SAMPLE="$( cat $1 )" # TS-9_leaf_Biseq_CHH
INDIR="/mnt/disk2/vibanez/descriptive-stat/aa_tmp-bed"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
#OUTDIR="/mnt/disk2/vibanez/descriptive-stat/ab_annotated-bed"
OUTDIR="/mnt/disk3/vibanez/descriptive-stat/ab_annotated-bed"
echo $SAMPLE
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t"  $INDIR/$SAMPLE.bed.tmp > $OUTDIR/$SAMPLE.bed

#echo "exon-TE"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/exon-TE.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_exon-TE.bed"

#echo "exon"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/exon.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_exon.bed"

#echo "intron-TE"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/intron-TE.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_intron-TE.bed"

#echo "intron"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/intron.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_intron.bed"

#echo "gene"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/gene.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_gene.bed"

#echo "gene-TE"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/gene-TE.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_gene-TE.bed"

#echo "TE-wo-gene"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/TE-wo-gene.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_TE-wo-gene.bed"

#echo "promotor"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/promotor.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_promotor.bed"

#echo "intergenic"
#intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $ANNODIR/intergenic.anno |\
#gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_intergenic.bed"
