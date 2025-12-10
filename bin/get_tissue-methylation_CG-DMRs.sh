SAMPLE="$( cat $1 )"
echo $SAMPLE
CHR="$( cat $1 | cut -d'_' -f5 )"
echo $CHR
SAMPLEPATH="/mnt/disk2/vibanez/01_bismark_analysis/04_processingFiles/ac_filter/bb_output"
LOCI_DIR="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/loci"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ad_merge_CG-DMR/aa_data"

#extrar la metilacion en cada ventana definida con todas las muestras
cat $SAMPLEPATH/$SAMPLE.filtered.bed | awk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
sortBed -i - | intersectBed -wb -a - -b $LOCI_DIR/"CG-DMR-loci_collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g' > $OUTPATH/$SAMPLE.methylation.bed

