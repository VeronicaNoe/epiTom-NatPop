CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
CG="/mnt/disk2/vibanez/02_methylkit/FRUITs/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/FRUITs/ah_merge_C-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ba_annotation/ab_out"

### C
echo '======= C-DMR ======'
echo '------500 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v > $OUTPATH/$CHR"_500-arribaTE.C-DMR.methylation"

echo '------ 500 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v > $OUTPATH/$CHR"_500-downTE.C-DMR.methylation"

echo '------1000 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v > $OUTPATH/$CHR"_1000-arribaTE.C-DMR.methylation"

echo '------1000 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v|\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v > $OUTPATH/$CHR"_1000-downTE.C-DMR.methylation"

echo '------1500 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v > $OUTPATH/$CHR"_1500-arribaTE.C-DMR.methylation"

echo '------1500 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v > $OUTPATH/$CHR"_1500-downTE.C-DMR.methylation"

echo '------2500 bp arriba TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"2500-arribaTE.anno"  |\
intersectBed -a - -b $ANNOPATH/"2500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"2500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno" -v > $OUTPATH/$CHR"_2500-arribaTE.C-DMR.methylation"


echo '------2500 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"2500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"2500-downTE.anno" |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"2500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno" -v > $OUTPATH/$CHR"_2500-downTE.C-DMR.methylation"


echo '======= C-DMR ======'
echo '------500 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v > $OUTPATH/$CHR"_500-arribaGENE.C-DMR.methylation"

echo '------ 500 bp down GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno"  > $OUTPATH/$CHR"_500-downGENE.C-DMR.methylation"

echo '------1000 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v > $OUTPATH/$CHR"_1000-arribaGENE.C-DMR.methylation"

echo '------1000 bp down TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v|\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno"  > $OUTPATH/$CHR"_1000-downGENE.C-DMR.methylation"

echo '------1500 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  > $OUTPATH/$CHR"_1500-arribaGENE.C-DMR.methylation"

echo '------1500 bp down GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno" |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v > $OUTPATH/$CHR"_1500-downGENE.C-DMR.methylation"

echo '------2500 bp arriba GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"2500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"2500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"2500-arribaGene.anno"  |\
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno" -v > $OUTPATH/$CHR"_2500-arribaGENE.C-DMR.methylation"

echo '------2500 bp down GENE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"2500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"2500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1500-downGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"2500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno"  > $OUTPATH/$CHR"_2500-downGENE.C-DMR.methylation"
