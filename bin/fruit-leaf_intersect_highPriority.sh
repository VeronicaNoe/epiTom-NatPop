TYPE="$( cat $1 | cut -d'_' -f2 )"
DMR="$( cat $1 | cut -d'_' -f1 )"
echo $TYPE

INDIR="/mnt/disk2/vibanez/02_methylkit/FRUITs/cb_intersections"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"

echo '------gene'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed"  -b $ANNOPATH/"gene.anno" -wa > $INDIR/$DMR.$TYPE."gene.methylation.bed"

echo '------gene plus TE'
intersectBed -a $INDIR/$DMR.$TYPE".loci.bed"  -b $ANNOPATH/"gene-TE.anno" -wa > $INDIR/$DMR.$TYPE."gene-TE.methylation.bed"

echo '------TE'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"gene.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno"  -wa > $INDIR/$DMR.$TYPE."TE-wo-gene.methylation.bed"

echo '------promotor'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"gene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"promotor.anno"  -wa > $INDIR/$DMR.$TYPE."promotor.methylation.bed"

echo '------intergenic'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"gene.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -v -wa | intersectBed -a - -b $ANNOPATH/"promotor.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"intergenic.anno" -wa > $INDIR/$DMR.$TYPE."intergenic.methylation.bed"

echo '------5UTR'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"5UTR.anno" -wa > $INDIR/$DMR.$TYPE."5UTR.methylation.bed"
echo '------3UTR'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"3UTR.anno" -wa > $INDIR/$DMR.$TYPE."3UTR.methylation.bed"
echo '------exon'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"promotor.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -wa > $INDIR/$DMR.$TYPE."exon.methylation.bed"
echo '------exon-TE'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"promotor.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon-TE.anno" -wa > $INDIR/$DMR.$TYPE."exon-TE.methylation.bed"
echo '------intron'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"promotor.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"intron.anno" -wa > $INDIR/$DMR.$TYPE."intron.methylation.bed"
echo '------intron-TE'
intersectBed -a $INDIR/$DMR.$TYPE."loci.bed" -b $ANNOPATH/"promotor.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon-TE.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"intron-TE.anno" -wa > $INDIR/$DMR.$TYPE."intron-TE.methylation.bed"
