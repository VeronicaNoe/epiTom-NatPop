CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
CG="/mnt/disk2/vibanez/02_methylkit/FRUITs/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/FRUITs/ah_merge_C-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ba_annotation/ab_out"
### CG
echo '======= CG-DMR ======'
echo '------gene'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -wa | bedtools sort -i - | mergeBed -i -  > $OUTPATH/$CHR"_gene.CG-DMR.methylation.tmp"
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_gene.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_gene.CG-DMR.methylation"

echo '------gene plus TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -wa | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_gene-TE.CG-DMR.methylation.tmp"
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_gene-TE.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_gene-TE.CG-DMR.methylation"

echo '------TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v  -wa  | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v  -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno"  -wa | bedtools sort -i - |\
mergeBed -i -  > $OUTPATH/$CHR"_TE-wo-gene.CG-DMR.methylation.tmp"

cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_TE-wo-gene.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_TE-wo-gene.CG-DMR.methylation"

echo '------promotor'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v   | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v |  mergeBed -i -  |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -v  | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"promotor.anno"  -wa | bedtools sort -i - |\
mergeBed -i - > $OUTPATH/$CHR"_promotor.CG-DMR.methylation.tmp"

cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_promotor.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_promotor.CG-DMR.methylation"

echo '------intergenic'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -v -wa| intersectBed -a - -b $ANNOPATH/"promotor.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"intergenic.anno" -wa | bedtools sort -i - > $OUTPATH/$CHR"_intergenic.CG-DMR.methylation.tmp"

cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_intergenic.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_intergenic.CG-DMR.methylation"

echo '------5UTR'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"5UTR-wo-TE.anno"  -wa > $OUTPATH/$CHR"_5UTR-wo-TE.CG-DMR.methylation"
echo '------3UTR'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"3UTR-wo-TE.anno"  -wa > $OUTPATH/$CHR"_3UTR-wo-TE.CG-DMR.methylation"
echo '------exon'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v  -wa | intersectBed -a - -b $ANNOPATH/"exon-wo-TE.anno" -wa> $OUTPATH/$CHR"_exon-wo-TE.CG-DMR.methylation"
echo '------intron'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v -wa | intersectBed -a - -b $ANNOPATH/"exon-wo-TE.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"intron-wo-TE.anno" -wa > $OUTPATH/$CHR"_intron-wo-TE.CG-DMR.methylation"

### C
echo '======= C-DMR ======'
echo '------gene'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -wa | bedtools sort -i - | mergeBed -i -  > $OUTPATH/$CHR"_gene.C-DMR.methylation.tmp"
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_gene.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_gene.C-DMR.methylation"

echo '------gene plus TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -wa | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_gene-TE.C-DMR.methylation.tmp"
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_gene-TE.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_gene-TE.C-DMR.methylation"

echo '------TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v  -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v  -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -wa | mergeBed -i -  > $OUTPATH/$CHR"_TE-wo-gene.C-DMR.methylation.tmp"

cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b  $OUTPATH/$CHR"_TE-wo-gene.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_TE-wo-gene.C-DMR.methylation"

echo '------promotor'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v  -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno" -v  -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"promotor.anno"  -wa | mergeBed -i - > $OUTPATH/$CHR"_promotor.C-DMR.methylation.tmp"

cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_promotor.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_promotor.C-DMR.methylation"


echo '------intergenic'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -v  -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -v -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.anno"  -v -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v -wa | mergeBed -i - |\
intersectBed -a - -b $ANNOPATH/"intergenic.anno" -wa | mergeBed -i - > $OUTPATH/$CHR"_intergenic.C-DMR.methylation.tmp"
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $OUTPATH/$CHR"_intergenic.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_intergenic.C-DMR.methylation"


echo '------5UTR'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"5UTR-wo-TE.anno"  -wa > $OUTPATH/$CHR"_5UTR-wo-TE.C-DMR.methylation"
echo '------3UTR'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"3UTR-wo-TE.anno"  -wa > $OUTPATH/$CHR"_3UTR-wo-TE.C-DMR.methylation"
echo '------exon'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v -wa | intersectBed -a - -b $ANNOPATH/"exon-wo-TE.anno" -wa > $OUTPATH/$CHR"_exon-wo-TE.C-DMR.methylation"
echo '------intron'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v -wa | intersectBed -a - -b $ANNOPATH/"exon-wo-TE.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"intron-wo-TE.anno"  -wa> $OUTPATH/$CHR"_intron-wo-TE.C-DMR.methylation"
