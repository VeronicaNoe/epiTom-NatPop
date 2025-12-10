CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
CG="/mnt/disk2/vibanez/02_methylkit/FRUITs/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/FRUITs/ah_merge_C-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ba_annotation/ab_out"

### CG
echo '======= CG-DMR ======'
echo '------500 bp arriba TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v > $OUTPATH/$CHR"_500-arribaTE.CG-DMR.methylation"

echo '------ 500 bp down TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v > $OUTPATH/$CHR"_500-downTE.CG-DMR.methylation"

echo '------1000 bp arriba TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v > $OUTPATH/$CHR"_1000-arribaTE.CG-DMR.methylation"

echo '------1000 bp down TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v|\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v > $OUTPATH/$CHR"_1000-downTE.CG-DMR.methylation"

echo '------1500 bp arriba TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v > $OUTPATH/$CHR"_1500-arribaTE.CG-DMR.methylation"

echo '------1500 bp down TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v > $OUTPATH/$CHR"_1500-downTE.CG-DMR.methylation"

echo '------2500 bp arriba TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno" -v > $OUTPATH/$CHR"_2500-arribaTE.CG-DMR.methylation"


echo '------2500 bp down TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno" -v > $OUTPATH/$CHR"_2500-downTE.CG-DMR.methylation"


echo '======= CG-DMR ======'
echo '------500 bp arriba GENE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v > $OUTPATH/$CHR"_500-arribaGENE.CG-DMR.methylation"

echo '------ 500 bp down GENE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno"  > $OUTPATH/$CHR"_500-downGENE.CG-DMR.methylation"

echo '------1000 bp arriba GENE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno" -v > $OUTPATH/$CHR"_1000-arribaGENE.CG-DMR.methylation"

echo '------1000 bp down TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"1000-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaTE.anno"  -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaTE.anno" -v|\
intersectBed -a - -b $ANNOPATH/"500-downTE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"500-downGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-arribaGene.anno" -v |\
intersectBed -a - -b $ANNOPATH/"1000-downGene.anno"  > $OUTPATH/$CHR"_1000-downGENE.CG-DMR.methylation"

echo '------1500 bp arriba GENE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  > $OUTPATH/$CHR"_1500-arribaGENE.CG-DMR.methylation"

echo '------1500 bp down GENE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"1500-arribaGene.anno"  -v > $OUTPATH/$CHR"_1500-downGENE.CG-DMR.methylation"

echo '------2500 bp arriba GENE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno" -v > $OUTPATH/$CHR"_2500-arribaGENE.CG-DMR.methylation"

echo '------2500 bp down GENE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | awk '{OFS="\t"}{if (NR>1){print} }' |\
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
intersectBed -a - -b $ANNOPATH/"2500-downGene.anno"  > $OUTPATH/$CHR"_2500-downGENE.CG-DMR.methylation"
