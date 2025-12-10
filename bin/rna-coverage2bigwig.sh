#!/bin/bash
SAMPLE="$( cat $1 )"
#INDIR="/mnt/data6/vibanez/IGV/data/RNAseq/KOs"
INDIR="/mnt/disk2/vibanez/07_rnaseq-processing/07.2_sorting"
#NAME="$( echo $SAMPLE | cut -d'.' -f1 | cut -d'_' -f1,2 )"

samtools sort -@ 20 -o ${INDIR}/${SAMPLE}_sorted.bam ${INDIR}/${SAMPLE}_merged.bam
echo "========= coverage ========="
samtools depth $INDIR/$SAMPLE.sorted.bam | grep 'SL2.50'| awk '{OFS="\t"}{print $1, $2,$2+1,$3}' > $INDIR/$SAMPLE.coverage.bed
echo "========= bigwig   ========="
# Sort bedGraph file
#sort -k1,1 -k2,2n $INDIR/$NAME.coverage.bed > $INDIR/$NAME.sorted.coverage.bed
~/bin/bedGraphToBigWig $INDIR/$NAME.coverage.bed ~/bin/chr.size.bed $INDIR/$NAME.bw



