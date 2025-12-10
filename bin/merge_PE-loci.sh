ANNO="$( cat $1 | cut -d'_' -f2)"
CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR "-" $ANNO

INDIR="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/ab_merged/tmp"
OUTDIR="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/ab_merged"
#ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/metaplot" # geneParts version
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data" # geneParts version
CHR_SIZE="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/chr.size.bed"
EPI="$( ls $INDIR/*$CHR"_"$ANNO.epiallele | tr '\n' '\t' )"
METH="$( ls $INDIR/*$CHR"_"$ANNO.methLevels | tr '\n' '\t' )" # geneParts version

#ls $INDIR/*$CHR"_"$ANNO.epiallele | sed 's./mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/ab_merged/tmp/..g' |\
#cut -d'_' -f1 > $OUTDIR/00_colNames.tsv

bedtools unionbedg -i $EPI -empty -g $CHR_SIZE -filler UM |\
intersectBed -a - -b $INDIR/$CHR"_"$ANNO".collapsed" -wb |\
intersectBed -a - -b $ANNODIR/allGenes.bed -wb > $OUTDIR/$CHR"_"$ANNO".epiallele"

#cat $OUTDIR/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }'| awk '{OFS="\t"}{print "chr","start","end", $0}'|\
#cat - $OUTDIR/$CHR"_"$ANNO".epiallele.tmp" |  sed 's/ /\t/g' > $OUTDIR/$CHR"_"$ANNO".epiallele"

#ls $METH
bedtools unionbedg -i $METH -empty -g $CHR_SIZE -filler UM |\
intersectBed -a - -b $INDIR/$CHR"_"$ANNO".collapsed" -wb |\
intersectBed -a - -b $ANNODIR/allGenes.bed -wb > $OUTDIR/$CHR"_"$ANNO".methLevels"

#cat $OUTDIR/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }'| awk '{OFS="\t"}{print "chr","start","end", $0}'|\
#cat - $OUTDIR/$CHR"_"$ANNO".methLevels.tmp" |  sed 's/ /\t/g' > $OUTDIR/$CHR"_"$ANNO".methLevels"
