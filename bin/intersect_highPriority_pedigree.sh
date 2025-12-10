SAMPLE="$( cat $1 )"
SAMPLE_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation"
ANNOPATH="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUTPATH="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection/bb_pedigree"

echo '------gene'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -wa | bedtools sort -i - | mergeBed -i - > $OUTPATH/${SAMPLE}"_gene.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$SAMPLE"_gene.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_gene.methylation"

echo '------gene plus TE'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -wa -v | intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -wa |\
bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_gene-TE.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_gene-TE.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_gene-TE.methylation"

echo '------TE'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v  -wa  | intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v  -wa |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno"  -wa | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_TE-wo-gene.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$SAMPLE"_TE-wo-gene.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_TE-wo-gene.methylation"

echo '------promotor'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v  | intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_promotor.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_promotor.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_promotor.methylation"

echo '------intergenic'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v -wa | intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v -wa  |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -v -wa| intersectBed -a - -b $ANNOPATH/"promotor.anno" -v  -wa|\
bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_intergenic.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_intergenic.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_intergenic.methylation"

echo '------5UTR'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno"  -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_5UTR.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_5UTR.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_5UTR.methylation"

echo '------3UTR'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"3UTR.anno"  -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_3UTR.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_3UTR.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_3UTR.methylation"

echo '------exon'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v  -wa |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_exon.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_exon.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_exon.methylation"

echo '------intron'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v -wa | intersectBed -a - -b $ANNOPATH/"exon.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"intron.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_intron.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_intron.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_intron.methylation"

echo '------exon-TE'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v |\
intersectBed -a - -b $ANNOPAHT/"exon-TE.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_exon-TE.CG-DMR.methylation.tmp"
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_exon-TE.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_exon-TE.methylation"

echo '------intron-TE'
cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v | intersectBed -a - -b $ANNOPATH/"exon.anno" -v |\
intersectBed -a - -b $ANNOPATH/"intron.anno" -v |
intersectBed -a - -b $ANNOPATH/"intron-TE.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$SAMPLE"_intron-TE.CG-DMR.methylation.tmp"

cat ${SAMPLE_DIR}/${SAMPLE}".merged.methylation.bed" | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$SAMPLE"_intron-TE.CG-DMR.methylation.tmp" > $OUTPATH/$SAMPLE"_intron-TE.methylation"


rm $OUTPATH/$SAMPLE*.tmp
