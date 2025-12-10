#!/bin/bash
CHR="$( cat $1 )"
echo $CHR
WD="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes"
INPUT=$WD"/aa_data"
OUTPUT=$WD"/bd_output"

INMETH="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/methylation/aa_data"
OUTMETH="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/methylation/bb_output"

INTRON_LIST="$( ls $INPUT/$CHR*"_epialleles-over-intron.bed" | grep -v 'TS-253' | tr '\n' '\t' )"
EXON_LIST="$( ls $INPUT/$CHR*"_epialleles-over-exon.bed" | grep -v 'TS-253' | tr '\n' '\t' )"
INTRON_TE_LIST="$( ls $INPUT/$CHR*"_epialleles-over-intron-TE.bed" | grep -v 'TS-253' | tr '\n' '\t' )"
EXON_TE_LIST="$( ls $INPUT/$CHR*"_epialleles-over-exon-TE.bed" | grep -v 'TS-253' | tr '\n' '\t' )"

INTRON_LIST_METH="$( ls $INMETH/$CHR*"_epialleles-over-intron.bed" | grep -v 'TS-253' |tr '\n' '\t' )"
EXON_LIST_METH="$( ls $INMETH/$CHR*"_epialleles-over-exon.bed" | grep -v 'TS-253' |tr '\n' '\t' )"
INTRON_TE_LIST_METH="$( ls $INMETH/$CHR*"_epialleles-over-intron-TE.bed" | grep -v 'TS-253' | tr '\n' '\t' )"
EXON_TE_LIST_METH="$( ls $INMETH/$CHR*"_epialleles-over-exon-TE.bed" | grep -v 'TS-253' | tr '\n' '\t' )"

CHR_SIZE="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ah_merge_C-DMR/chr.size.bed"
echo "Merging epialleles"
bedtools unionbedg -i $INTRON_LIST -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/intron-wo-TE.merged.anno  > $OUTPUT/$CHR"_intron.tmp"
bedtools unionbedg -i $EXON_LIST -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/exon-wo-TE.merged.anno  > $OUTPUT/$CHR"_exon.tmp"
bedtools unionbedg -i $INTRON_TE_LIST -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/intron-TE.merged.anno  > $OUTPUT/$CHR"_intronTE.tmp"
bedtools unionbedg -i $EXON_TE_LIST -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/exon-TE.merged.anno  > $OUTPUT/$CHR"_exonTE.tmp"
echo "Merging methylation"
bedtools unionbedg -i $INTRON_LIST_METH -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/intron-wo-TE.merged.anno  > $OUTMETH/$CHR"_intron.tmp"
bedtools unionbedg -i $EXON_LIST_METH -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/exon-wo-TE.merged.anno  > $OUTMETH/$CHR"_exon.tmp"
bedtools unionbedg -i $INTRON_TE_LIST_METH -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/intron-TE.merged.anno  > $OUTMETH/$CHR"_intronTE.tmp"
bedtools unionbedg -i $EXON_TE_LIST_METH -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b /mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data/exon-TE.merged.anno  > $OUTMETH/$CHR"_exonTE.tmp"
