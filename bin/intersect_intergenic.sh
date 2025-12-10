CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ah_merge_C-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/ab_out"

### CG
echo '======= CG-DMR ======'
echo '------intergenic'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed" | awk '{OFS="\t"}{if (NR>1){print} }' |intersectBed -a - -b $ANNOPATH/"intergenic.anno" -wa | bedtools sort -i - > $OUTPATH/$CHR"_intergenic.CG-DMR.methylation"
### C
echo '======= C-DMR ======'
echo '------intergenic'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' | intersectBed -a - -b $ANNOPATH/"intergenic.anno" -wa | bedtools sort -i -  > $OUTPATH/$CHR"_intergenic.C-DMR.methylation"
