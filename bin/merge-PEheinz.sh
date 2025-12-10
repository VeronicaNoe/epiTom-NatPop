#!/bin/bash
WD="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes"
INPUT=$WD"/aa_data"
OUTPUT=$WD"/bd_output"

INMETH="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/methylation/aa_data"
OUTMETH="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/methylation/bb_output"

ls $INPUT/ch01*"_epialleles-over-exon.bed" | grep -v 'TS-253' |\
sed 's./mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/aa_data/..g' |\
sed 's/_epialleles-over-exon.bed//g' | cut -d'_' -f2  > $OUTPUT/00_colNames.tsv

ls $INMETH/ch01*"_epialleles-over-exon.bed" | grep -v 'TS-253' |\
sed 's./mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/methylation/aa_data/..g' |\
sed 's/_epialleles-over-exon.bed//g' | cut -d'_' -f2  > $OUTMETH/00_colNames.tsv

# add names
cat $OUTPUT/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }' | awk '{ OFS="\t" } { print "chr", "start", "end", $0}' | cat - $OUTPUT/*"exon.tmp" > $OUTPUT/PEexon.bed
cat $OUTPUT/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }' | awk '{ OFS="\t" } { print "chr", "start", "end", $0}' | cat - $OUTPUT/*"intron.tmp" > $OUTPUT/PEintron.bed
cat $OUTPUT/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }' | awk '{ OFS="\t" } { print "chr", "start", "end", $0}' | cat - $OUTPUT/*"intronTE.tmp" > $OUTPUT/PEintronTE.bed
cat $OUTPUT/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }' | awk '{ OFS="\t" } { print "chr", "start", "end", $0}' | cat - $OUTPUT/*"_exonTE.tmp" > $OUTPUT/PEexonTE.bed

echo "TS-253"
cat $OUTPUT/PEintron.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INPUT/"TS-253_epialleles-over-intron.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTPUT/"TS-253_intron.tmp"
cat $OUTPUT/PEexon.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INPUT/"TS-253_epialleles-over-exon.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTPUT/"TS-253_exon.tmp"
cat $OUTPUT/PEintronTE.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INPUT/"TS-253_epialleles-over-intron-TE.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTPUT/"TS-253_intronTE.tmp"
cat $OUTPUT/PEexonTE.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INPUT/"TS-253_epialleles-over-exon-TE.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTPUT/"TS-253_exonTE.tmp"
cat $OUTPUT/PEintron.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INMETH/"TS-253_epialleles-over-intron.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTMETH/"TS-253_intron.tmp"
cat $OUTPUT/PEexon.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INMETH/"TS-253_epialleles-over-exon.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTMETH/"TS-253_exon.tmp"
cat $OUTPUT/PEintronTE.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INMETH/"TS-253_epialleles-over-intron-TE.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTMETH/"TS-253_intronTE.tmp"
cat $OUTPUT/PEexonTE.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | mergeBed -i - | intersectBed -a - -b $INMETH/"TS-253_epialleles-over-exon-TE.bed" -wao |\
awk '{OFS="\t"} {print $1, $2, $3, $7}' > $OUTMETH/"TS-253_exonTE.tmp"

cat $OUTPUT/00_TS-253.tsv $OUTPUT/TS-253_exon.tmp > $OUTPUT/TS-253_exon.bed
cat $OUTPUT/00_TS-253.tsv $OUTPUT/TS-253_exonTE.tmp > $OUTPUT/TS-253_exonTE.bed
cat $OUTPUT/00_TS-253.tsv $OUTPUT/TS-253_intron.tmp > $OUTPUT/TS-253_intron.bed
cat $OUTPUT/00_TS-253.tsv $OUTPUT/TS-253_intronTE.tmp > $OUTPUT/TS-253_intronTE.bed

cat $OUTPUT/PEexon.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/exonID.merged -wao |awk '{OFS="\t"}{ print $1,$2,$3, $7}'  > IDexon.bed
cat $OUTPUT/PEintron.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - |intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/intronID.merged -wao |awk '{OFS="\t"}{ print $1,$2,$3, $7}' > IDintron.bed
cat $OUTPUT/PEexonTE.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/exonID.merged -wao |awk '{OFS="\t"}{ print $1,$2,$3, $7}'  > IDexonTE.bed
cat $OUTPUT/PEintronTE.bed | awk '{OFS="\t"} NR>1 {print $1, $2, $3}' | sortBed -i - |intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/intronID.merged -wao | awk '{OFS="\t"}{ print $1,$2,$3, $7}' > IDintronTE.bed
