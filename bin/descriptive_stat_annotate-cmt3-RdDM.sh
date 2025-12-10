#!/bin/bash
echo 'Activate gawk environment'
CTXT="$( cat $1 )"
INDIR="/mnt/disk2/vibanez/descriptive-stat/aa_tmp-bed"
KODIR="/mnt/disk2/vibanez/04_natural-experimental_DMR/aa_KO"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTDIR="/mnt/disk2/vibanez/04_natural-experimental_DMR/annotate_kos"
intersectBed -a $KODIR/"cmt3-kyp_"$CTXT"-DMR.methylation.bed" -b $KODIR/"slnrp-de_"$CTXT"-DMR.methylation.bed" > $KODIR/"cmt3-RdDM"_$CTXT".bed"

echo "exon-TE"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/exon-TE.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_exon-TE.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/exon-TE.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_exon-TE.bed"

echo "exon"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/exon.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_exon.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/exon.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_exon.bed"

echo "intron-TE"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/intron-TE.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_intron-TE.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/intron-TE.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_intron-TE.bed"

echo "intron"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/intron.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_intron.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/intron.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_intron.bed"

echo "gene"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/gene.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_gene.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/gene.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_gene.bed"

echo "gene-TE"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/gene-TE.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_gene-TE.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/gene-TE.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_gene-TE.bed"

echo "TE-wo-gene"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/TE-wo-gene.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_TE-wo-gene.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/TE-wo-gene.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_TE-wo-gene.bed"

echo "promotor"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/promotor.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_promotor.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/promotor.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_promotor.bed"


echo "intergenic"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/intergenic.anno |\
gawk '{sum5 += $5; sum6 += $6} END {print sum5, sum6, sum6/sum5}' OFS="\t" > $OUTDIR/"control_cmt3-RdDM"_$CTXT"_intergenic.bed"
intersectBed -a $KODIR/"cmt3-RdDM"_$CTXT".bed" -b $ANNODIR/intergenic.anno |\
gawk '{sum8 += $8; sum9 += $9} END {print sum8, sum9, sum9/sum8}' OFS="\t" > $OUTDIR/"KO_cmt3-RdDM"_$CTXT"_intergenic.bed"
