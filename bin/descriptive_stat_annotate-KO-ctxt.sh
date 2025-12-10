#!/bin/bash
echo 'Activate gawk environment'
SAMPLE="$( cat $1 | cut -d'_' -f1,2,3,4)" # TS-9_leaf_Biseq_CHH
KO="$( cat $1 | cut -d'_' -f5,6)"
INDIR="/mnt/disk2/vibanez/descriptive-stat/aa_tmp-bed"
KODIR="/mnt/disk2/vibanez/04_natural-experimental_DMR/aa_KO"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
#OUTDIR="/mnt/disk2/vibanez/descriptive-stat/ab_annotated-bed"
OUTDIR="/mnt/disk3/vibanez/descriptive-stat/ab_annotated-bed"
echo $SAMPLE "-" $KO
intersectBed -a $INDIR/$SAMPLE.bed.tmp -b $KODIR/$KO".ctxt-merged" -wb | gawk '{OFS="\t"}{ if ($18<0) {print}}' > $OUTDIR/$SAMPLE"_"$KO".methylation"
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t"  $OUTDIR/$SAMPLE"_"$KO".methylation" > $OUTDIR/$SAMPLE"_"$KO".bed"

echo "exon-TE"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/exon-TE.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_exon-TE.bed"

echo "exon"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/exon.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_exon.bed"

echo "intron-TE"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/intron-TE.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_intron-TE.bed"

echo "intron"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/intron.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_intron.bed"

echo "gene"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/gene.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_gene.bed"

echo "gene-TE"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/gene-TE.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_gene-TE.bed"

echo "TE-wo-gene"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/TE-wo-gene.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_TE-wo-gene.bed"

echo "promotor"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/promotor.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_promotor.bed"

echo "intergenic"
intersectBed -a $OUTDIR/$SAMPLE"_"$KO".methylation" -b $ANNODIR/intergenic.anno |\
gawk '{sum4 += $4; sum5 += $5} END {print sum4, sum5, sum4/(sum4+sum5)}' OFS="\t" > $OUTDIR/$SAMPLE"_"$KO"_intergenic.bed"
