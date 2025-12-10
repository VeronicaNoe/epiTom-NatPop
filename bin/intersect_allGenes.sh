CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
DMRs="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ac_merge-methylation/bb_output"
ANNO_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUT_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection"
### CG
echo '======= CG-DMR ======'
echo '------gene'
cat $DMRs/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNO_DIR/"allGenes.bed" -wa | bedtools sort -i - | mergeBed -i - > $OUT_DIR/$CHR"_genes.CG-DMR.methylation.tmp"
cat $DMRs/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUT_DIR/$CHR"_genes.CG-DMR.methylation.tmp" > $OUT_DIR/$CHR"_genes.CG-DMR.methylation"

### C
echo '======= C-DMR ======'
echo '------gene'
        cat $DMRs/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNO_DIR/"allGenes.bed" -wa | bedtools sort -i - | mergeBed -i - > $OUT_DIR/$CHR"_genes.C-DMR.methylation.tmp"
        cat $DMRs/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUT_DIR/$CHR"_genes.C-DMR.methylation.tmp" > $OUT_DIR/$CHR"_genes.C-DMR.methylation"

rm $OUT_DIR/$CHR*.tmp
