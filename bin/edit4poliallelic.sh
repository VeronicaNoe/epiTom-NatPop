#!/bin/bash
WD="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes"
OUTPUT="aa_data"

#intersectBed -a $1 -b $WD/genes-C-DMR.bed | awk '{OFS="\t"}{print $1, $2, $3, $4=1}' > $WD/$OUTPUT/$1
intersectBed -a $1 -b $WD/genes-PE.bed | awk '{OFS="\t"}{print $1, $2, $3, $4=1}' > $WD/$OUTPUT/$1
