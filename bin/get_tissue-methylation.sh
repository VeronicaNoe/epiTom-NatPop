SAMPLE="$( cat $1 | cut -d'_' -f1,2)"
echo $SAMPLE
CHR="$( cat $1 | cut -d'_' -f3 )"
echo $CHR
SAMPLEPATH="/mnt/disk2/vibanez/01_bismark_analysis/04_processingFiles/ac_filter/bb_output"
LOCI_DIR="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge-CG/aa_loci"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge-CG/ab_methBySample/bb_output"

#extrar la metilacion en cada ventana definida con todas las muestras
#CG-DMR
cat $SAMPLEPATH/$SAMPLE*"CG_"$CHR".filtered.bed" | awk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
intersectBed -wb -a - -b $LOCI_DIR/$CHR"_CG-DMR-loci_collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g' > $OUTPATH/$SAMPLE"_"$CHR"_CG-DMR.methylation.bed"

#C-DMR
cat $SAMPLEPATH/$SAMPLE*"CG_"$CHR".filtered.bed" | awk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
intersectBed -wb -a - -b $LOCI_DIR/$CHR"_C-DMR-loci_collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' > $OUTPATH/$SAMPLE"_"$CHR"_CG_C-DMR.methylation.tmp"

cat $SAMPLEPATH/$SAMPLE*"CHG_"$CHR".filtered.bed" | awk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
intersectBed -wb -a - -b $LOCI_DIR/$CHR"_C-DMR-loci_collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' > $OUTPATH/$SAMPLE"_"$CHR"_CHG_C-DMR.methylation.tmp"

cat $SAMPLEPATH/$SAMPLE*"CHH_"$CHR".filtered.bed" | awk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
intersectBed -wb -a - -b $LOCI_DIR/$CHR"_C-DMR-loci_collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' > $OUTPATH/$SAMPLE"_"$CHR"_CHH_C-DMR.methylation.tmp"

# get weigth methylation C-DMR
cat $OUTPATH/$SAMPLE"_"$CHR*"_C-DMR.methylation.tmp" | sort -k 1.4n -k 2n |awk '{OFS="\t"}{print $1,$2,$3,$4,$5}' |\
mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g'> $OUTPATH/$SAMPLE"_"$CHR"_C-DMR.methylation.bed"

rm $OUTPATH/$SAMPLE"_"$CHR*"_C-DMR.methylation.tmp"
