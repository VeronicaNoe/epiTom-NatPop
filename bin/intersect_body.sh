CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ah_merge_C-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/gbM"
### CG
echo '======= CG-DMR ======'
echo '------q1 gene'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q1-gene.anno" > $OUTPATH/$CHR"_q1-gene.CG-DMR.methylation"
echo '------q2 gene'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q2-gene.anno" > $OUTPATH/$CHR"_q2-gene.CG-DMR.methylation"
echo '------q3 gene'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q3-gene.anno" > $OUTPATH/$CHR"_q3-gene.CG-DMR.methylation"
echo '------q4 gene'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q4-gene.anno" > $OUTPATH/$CHR"_q4-gene.CG-DMR.methylation"

echo '------q1 TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q1-TE.anno" > $OUTPATH/$CHR"_q1-TE.CG-DMR.methylation"
echo '------q2 TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q2-TE.anno" > $OUTPATH/$CHR"_q2-TE.CG-DMR.methylation"
echo '------q3 TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q3-TE.anno" > $OUTPATH/$CHR"_q3-TE.CG-DMR.methylation"
echo '------q4 TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q4-TE.anno" > $OUTPATH/$CHR"_q4-TE.CG-DMR.methylation"


### C
echo '======= C-DMR ======'
echo '------q1 gene'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q1-gene.anno" > $OUTPATH/$CHR"_q1-gene.C-DMR.methylation"
echo '------q2 gene'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q2-gene.anno" > $OUTPATH/$CHR"_q2-gene.C-DMR.methylation"
echo '------q3 gene'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q3-gene.anno" > $OUTPATH/$CHR"_q3-gene.C-DMR.methylation"
echo '------q4 gene'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q4-gene.anno" > $OUTPATH/$CHR"_q4-gene.C-DMR.methylation"

echo '------q1 TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q1-TE.anno" > $OUTPATH/$CHR"_q1-TE.C-DMR.methylation"
echo '------q2 TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q2-TE.anno" > $OUTPATH/$CHR"_q2-TE.C-DMR.methylation"
echo '------q3 TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q3-TE.anno" > $OUTPATH/$CHR"_q3-TE.C-DMR.methylation"
echo '------q4 TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"q4-TE.anno" > $OUTPATH/$CHR"_q4-TE.C-DMR.methylation"
