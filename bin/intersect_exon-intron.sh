CHR="$( cat $1 | cut -d'_' -f1 )"
echo $CHR
CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ah_merge_C-DMR/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/ab_out"
### CG
echo '======= CG-DMR ======'
echo '------exon'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v  -wa |\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_exon.CG-DMR.methylation.tmp"
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_exon.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_exon.CG-DMR.methylation"

echo '------intron'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"intron.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_intron.CG-DMR.methylation.tmp"
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_intron.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_intron.CG-DMR.methylation"

echo '------exon-TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v |\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon-TE.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_exon-TE.CG-DMR.methylation.tmp"
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_exon-TE.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_exon-TE.CG-DMR.methylation"

echo '------intron-TE'
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -v |\
intersectBed -a - -b $ANNOPATH/"exon-TE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"intron.anno" -v |\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"intron-TE.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_intron-TE.CG-DMR.methylation.tmp"
cat $CG/$CHR"_CG-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_intron-TE.CG-DMR.methylation.tmp" > $OUTPATH/$CHR"_intron-TE.CG-DMR.methylation"


### C
echo '======= C-DMR ======'
echo '------exon'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v  -wa |\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_exon.C-DMR.methylation.tmp"
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_exon.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_exon.C-DMR.methylation"

echo '------intron'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -v  -wa|\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"intron.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_intron.C-DMR.methylation.tmp"
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_intron.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_intron.C-DMR.methylation"

echo '------exon-TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v |\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"exon-TE.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_exon-TE.C-DMR.methylation.tmp"
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_exon-TE.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_exon-TE.C-DMR.methylation"

echo '------intron-TE'
cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNOPATH/"promotor.anno" -v |\
intersectBed -a - -b $ANNOPATH/"exon.anno" -v |\
intersectBed -a - -b $ANNOPATH/"exon-TE.anno" -v |\
intersectBed -a - -b $ANNOPATH/"intron.anno" -v |\
intersectBed -a - -b $ANNOPATH/"3UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"5UTR.anno" -v -wa |\
intersectBed -a - -b $ANNOPATH/"intron-TE.anno" -wa  | bedtools sort -i - | mergeBed -i - > $OUTPATH/$CHR"_intron-TE.C-DMR.methylation.tmp"

cat $C/$CHR"_C-DMR.merged.methylation.bed"  | gawk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a -  -b $OUTPATH/$CHR"_intron-TE.C-DMR.methylation.tmp" > $OUTPATH/$CHR"_intron-TE.C-DMR.methylation"
