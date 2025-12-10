#!/bin/bash
CHR="$( cat $1 )"
echo $CHR
SAMPLEPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ac_loci_C-DMRs/bb_output"
TMPFILES="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ac_loci_C-DMRs/data-tmp"
LOCI="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/loci"
# creating a small data set
for file in $(ls $SAMPLEPATH/$CHR*.C-DMR.bed | sed 's./mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/aa_data/..g' ); do
	cat $SAMPLEPATH/$file | awk '{OFS="\t"}{print $1,$2,$3,$6/$5}' | sed 's/,/./g' > $TMPFILES/$file
done
#count the overlaps
cat $TMPFILES/$CHR*.C-DMR.bed | sort -k 1.4n -k 2n |awk '{OFS="\t"}{print $1,$2,$3,$4}'|\
mergeBed -c 4,4 -o count,distinct | awk '{OFS="\t"}{print $1,$2,$3,$4,$5}'| uniq >$LOCI/$CHR"_C-DMR-loci_collapsed.bed"

#loop over all the files to have exactly the same coordinates for all, needed to run multiinter
for FILE in $(ls $SAMPLEPATH/$CHR*.C-DMR.bed | sed 's./mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/aa_data/..g' ); do
	awk '{OFS="\t"}{print $1,$2,$3}' $SAMPLEPATH/$FILE | intersectBed -wa -a $LOCI/$CHR"_C-DMR-loci_collapsed.bed" -b - -nonamecheck | uniq > $LOCI/${FILE%%.*}.C-DMR.merged.bed
done
#build the matrix
bedtools multiinter -i $LOCI/$CHR"_C-DMR-loci_collapsed.bed" $LOCI/$CHR*.C-DMR.merged.bed | awk '{print $1,$2,$3,$4,$5}' |awk '{OFS="\t"}{$4=$4-1}1'| uniq > $OUTPATH/$CHR"_C-DMR-loci_matrix.tsv"
#rm $TMP/$CHR*.merged.bed
