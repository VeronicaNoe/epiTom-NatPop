CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
C="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ac_merge-methylation/bb_output"
ANNOPATH="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUTPATH="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection"

### C
echo '======= TEs ======'
echo '------500 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v > $OUTPATH/$CHR"_TE-wo-gene_arriba-0500.C-DMR.methylation"

echo '------ 500 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-0500.anno" |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v > $OUTPATH/$CHR"_TE-wo-gene_down-0500.C-DMR.methylation"

echo '------1000 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v > $OUTPATH/$CHR"_TE-wo-gene_arriba-1000.C-DMR.methylation"

echo '------1000 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v|\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v > $OUTPATH/$CHR"_TE-wo-gene_down-1000.C-DMR.methylation"

echo '------1500 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno"  -v > $OUTPATH/$CHR"_TE-wo-gene_arriba-1500.C-DMR.methylation"

echo '------1500 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno"  -v > $OUTPATH/$CHR"_TE-wo-gene_down-1500.C-DMR.methylation"

echo '------2500 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-2500.anno"  |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-2500.anno" -v > $OUTPATH/$CHR"_TE-wo-gene_arriba-2500.C-DMR.methylation"


echo '------2500 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-2500.anno" |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-2500.anno" -v > $OUTPATH/$CHR"_TE-wo-gene_down-2500.C-DMR.methylation"

echo '======= GENE ======'
echo '------500 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v > $OUTPATH/$CHR"_gene_arriba-0500.C-DMR.methylation"

echo '------ 500 bp down GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" > $OUTPATH/$CHR"_gene_down-0500.C-DMR.methylation"

echo '------1000 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v > $OUTPATH/$CHR"_gene_arriba-1000.C-DMR.methylation"

echo '------1000 bp down GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" > $OUTPATH/$CHR"_gene_down-1000.C-DMR.methylation"

echo '------1500 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno" > $OUTPATH/$CHR"_gene_arriba-1500.C-DMR.methylation"

echo '------1500 bp down GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno"  -v > $OUTPATH/$CHR"_gene_down-1500.C-DMR.methylation"

echo '------2500 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-2500.anno" |\
intersectBed -a - -b $ANNOPATH/"gene_down-2500.anno" -v > $OUTPATH/$CHR"_gene_arriba-2500.C-DMR.methylation"

echo '------2500 bp down GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-0500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1000.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-1500.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"gene_arriba-2500.anno" -v |\
intersectBed -a - -b $ANNOPATH/"gene_down-2500.anno" > $OUTPATH/$CHR"_gene_down-2500.C-DMR.methylation"
