CHR="chrNonRef"
echo $CHR
CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ah_merge_C-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/NonRefSeq/aa_annotation/ba_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/NonRefSeq/aa_annotation/bb_output"
### CG
echo '======= CG-DMR ======'
echo '------gene'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"nonRef-genes.anno" -wa > $OUTPATH/$CHR"_gene.CG-DMR.methylation"
echo '------5UTR'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/tmp/"nonRef-5UTR.anno" -wa > $OUTPATH/$CHR"_5UTR.CG-DMR.methylation"
echo '------3UTR'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/tmp/"nonRef-3UTR.anno" -wa > $OUTPATH/$CHR"_3UTR.CG-DMR.methylation"
echo '------exon'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"nonRef-exon.anno" -wa > $OUTPATH/$CHR"_exon.CG-DMR.methylation"

cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"nonRef-genes.anno" -v > $OUTPATH/$CHR"_nonGene.CG-DMR"


### C
echo '======= C-DMR ======'
echo '------gene'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"nonRef-genes.anno" -wa > $OUTPATH/$CHR"_gene.C-DMR.methylation"
echo '------5UTR'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/tmp/"nonRef-5UTR.anno" -wa > $OUTPATH/$CHR"_5UTR.C-DMR.methylation"
echo '------3UTR'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/tmp/"nonRef-3UTR.anno" -wa > $OUTPATH/$CHR"_3UTR.C-DMR.methylation"
echo '------exon'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"nonRef-exon.anno" -wa > $OUTPATH/$CHR"_exon.C-DMR.methylation"
