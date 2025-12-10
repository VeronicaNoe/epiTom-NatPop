ANNO="$( cat $1 | cut -d'_' -f2)"
CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR "-" $ANNO

INDIR="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/aa_data"
OUTDIR="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/ab_merged/tmp"
#count the overlaps by each windows
cat $INDIR/*"_"$CHR"_epialleles-over-"$ANNO.bed | sort -k 1.4n -k 2n | gawk '{OFS="\t"}{print $1,$2,$3,$5}'|\
mergeBed -i - -c 4,4 -o distinct,count > $OUTDIR/$CHR"_"$ANNO".collapsed"

for FILE in $(ls $INDIR/*"_"$CHR"_epialleles-over-"$ANNO.bed | sed 's./mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/aa_data/..g' ); do
      SAMPLE="$( echo $FILE | cut -d'_' -f1)"
        #echo $SAMPLE
	cat $INDIR/$SAMPLE"_"$CHR"_epialleles-over-"$ANNO.bed | awk '{OFS="\t"}{print $1,$2,$3,$4}' |\
	intersectBed -a $OUTDIR/$CHR"_"$ANNO".collapsed" -b - -nonamecheck -wb | awk '{OFS="\t"}{print $1,$2,$3,$9}'> $OUTDIR/$SAMPLE"_"$CHR"_"$ANNO.methLevels
	cat $INDIR/$SAMPLE"_"$CHR"_epialleles-over-"$ANNO.bed | awk '{OFS="\t"}{print $1,$2,$3,$5}' |\
	intersectBed -a $OUTDIR/$CHR"_"$ANNO".collapsed" -b - -nonamecheck -wb | awk '{OFS="\t"}{print $1,$2,$3,$9}'> $OUTDIR/$SAMPLE"_"$CHR"_"$ANNO.epiallele
done
