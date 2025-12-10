#!/bin/bash
WD="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/aa_data"
OUTPUT="$WD/tmp"
ANNO="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ac_CG-DMR/bb_output/Heinz"
C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ad_C-DMR/bb_output/Heinz"
OUTMETH="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/methylation/aa_data"
## exon
cat $CG/*"TS-253_leaf.CG-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"exon-wo-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum > $OUTPUT/"TS-253.CG-exon"
cat $C/*"TS-253_leaf.C-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"exon-wo-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum > $OUTPUT/"TS-253.C-exon"
cat $CG/*"TS-253_leaf.CG-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"exon-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum > $OUTPUT/"TS-253.CG-exon-TE"
cat $C/*"TS-253_leaf.C-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"exon-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum  > $OUTPUT/"TS-253.C-exon-TE"

intersectBed -a $OUTPUT/"TS-253.CG-exon" -b $OUTPUT/"TS-253.C-exon" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/"TS-253_PE-exon.bed"
intersectBed -a $OUTPUT/"TS-253.CG-exon" -b $OUTPUT/"TS-253_PE-exon.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/"TS-253_only-CG-exon.bed"
intersectBed -a $OUTPUT/"TS-253.C-exon" -b $OUTPUT/"TS-253_PE-exon.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/"TS-253_only-C-exon.bed"

intersectBed -a $OUTPUT/"TS-253.CG-exon-TE" -b $OUTPUT/"TS-253.C-exon-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/"TS-253_PE-exon-TE.bed"
intersectBed -a $OUTPUT/"TS-253.CG-exon-TE" -b $OUTPUT/"TS-253_PE-exon-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/"TS-253_only-CG-exon-TE.bed"
intersectBed -a $OUTPUT/"TS-253.C-exon-TE" -b $OUTPUT/"TS-253_PE-exon-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/"TS-253_only-C-exon-TE.bed"

intersectBed -a $OUTPUT/"TS-253.CG-exon" -b $OUTPUT/"TS-253.C-exon" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/"TS-253_PE-exon.bed"
intersectBed -a $OUTPUT/"TS-253.CG-exon" -b $OUTPUT/"TS-253_PE-exon.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-CG-exon.bed"
intersectBed -a $OUTPUT/"TS-253.C-exon" -b $OUTPUT/"TS-253_PE-exon.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-C-exon.bed"
intersectBed -a $OUTPUT/"TS-253.CG-exon-TE" -b $OUTPUT/"TS-253.C-exon-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/"TS-253_PE-exon-TE.bed"
intersectBed -a $OUTPUT/"TS-253.CG-exon-TE" -b $OUTPUT/"TS-253_PE-exon-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-CG-exon-TE.bed"
intersectBed -a $OUTPUT/"TS-253.C-exon-TE" -b $OUTPUT/"TS-253_PE-exon-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-C-exon-TE.bed"

#intron
cat $CG/*"TS-253_leaf.CG-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"intron-wo-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum > $OUTPUT/"TS-253.CG-intron"
cat $C/*"TS-253_leaf.C-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"intron-wo-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum > $OUTPUT/"TS-253.C-intron"
cat $CG/*"TS-253_leaf.CG-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"intron-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum > $OUTPUT/"TS-253.CG-intron-TE"
cat $C/*"TS-253_leaf.C-DMR.bed" | sortBed -i | awk '{OFS="\t"}{print $1,$2,$3, $6/$5}' | sed 's/,/./g' |\
intersectBed -b - -a $ANNO/"intron-TE.merged.anno" -wao | awk '{if ($8>=99) {print}}'| mergeBed -i - -c 7,8 -o mean,sum  > $OUTPUT/"TS-253.C-intron-TE"

intersectBed -a $OUTPUT/"TS-253.CG-intron" -b $OUTPUT/"TS-253.C-intron" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/"TS-253_PE-intron.bed"
intersectBed -a $OUTPUT/"TS-253.CG-intron" -b $OUTPUT/"TS-253_PE-intron.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/"TS-253_only-CG-intron.bed"
intersectBed -a $OUTPUT/"TS-253.C-intron" -b $OUTPUT/"TS-253_PE-intron.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/"TS-253_only-C-intron.bed"

intersectBed -a $OUTPUT/"TS-253.CG-intron-TE" -b $OUTPUT/"TS-253.C-intron-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/"TS-253_PE-intron-TE.bed"
intersectBed -a $OUTPUT/"TS-253.CG-intron-TE" -b $OUTPUT/"TS-253_PE-intron-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/"TS-253_only-CG-intron-TE.bed"
intersectBed -a $OUTPUT/"TS-253.C-intron-TE" -b $OUTPUT/"TS-253_PE-intron-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/"TS-253_only-C-intron-TE.bed"

intersectBed -a $OUTPUT/"TS-253.CG-intron" -b $OUTPUT/"TS-253.C-intron" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/"TS-253_PE-intron.bed"
intersectBed -a $OUTPUT/"TS-253.CG-intron" -b $OUTPUT/"TS-253_PE-intron.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-CG-intron.bed"
intersectBed -a $OUTPUT/"TS-253.C-intron" -b $OUTPUT/"TS-253_PE-intron.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-C-intron.bed"
intersectBed -a $OUTPUT/"TS-253.CG-intron-TE" -b $OUTPUT/"TS-253.C-intron-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/"TS-253_PE-intron-TE.bed"
intersectBed -a $OUTPUT/"TS-253.CG-intron-TE" -b $OUTPUT/"TS-253_PE-intron-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-CG-intron-TE.bed"
intersectBed -a $OUTPUT/"TS-253.C-intron-TE" -b $OUTPUT/"TS-253_PE-intron-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/"TS-253_only-C-intron-TE.bed"

#
cat $OUTPUT/"TS-253"*"-exon.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $WD/"TS-253_epialleles-over-exon.bed"
cat $OUTPUT/"TS-253"*"-exon-TE.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $WD/"TS-253_epialleles-over-exon-TE.bed"
cat $OUTMETH/"TS-253"*"-exon.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $OUTMETH/"TS-253_epialleles-over-exon.bed"
cat $OUTMETH/"TS-253"*"-exon-TE.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $OUTMETH/"TS-253_epialleles-over-exon-TE.bed"

cat $OUTPUT/"TS-253"*"-intron.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $WD/"TS-253_epialleles-over-intron.bed"
cat $OUTPUT/"TS-253"*"-intron-TE.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $WD/"TS-253_epialleles-over-intron-TE.bed"
cat $OUTMETH/"TS-253"*"-intron.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $OUTMETH/"TS-253_epialleles-over-intron.bed"
cat $OUTMETH/"TS-253"*"-intron-TE.bed" | sortBed -i - |awk '{OFS="\t"}{print}' > $OUTMETH/"TS-253_epialleles-over-intron-TE.bed"
