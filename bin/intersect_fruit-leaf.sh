CG_LEAF="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ae_loci_CG-DMRs/tmp"
C_LEAF="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/af_loci_C-DMRs/tmp"
CG_FRUIT="/mnt/disk2/vibanez/02_methylkit/FRUITs/ad_loci_CG-DMRs/tmp"
C_FRUIT="/mnt/disk2/vibanez/02_methylkit/FRUITs/ae_loci_C-DMRs/tmp"
MERGED="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/loci"
CG_METH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ad_merge_CG-DMR/bb_output"
C_METH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ae_merge_C-DMR/bb_output"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/cb_intersections/"

#echo '======= CG-DMR ======'
#intersectBed -a $CG_LEAF/"CG-DMR-loci_collapsed.bed" -b $CG_FRUIT/"CG-DMR-loci_collapsed.bed" |\
#awk '{OFS="\t"} {print $1, $2, $3}' > $OUTPATH/"CG-DMR.common.loci"

#intersectBed -a $MERGED/"CG-DMR-loci_collapsed.bed" -b $CG_FRUIT/"CG-DMR-loci_collapsed.bed" -v |\
#awk '{OFS="\t"} {print $1, $2, $3}' > $OUTPATH/"CG-DMR.onlyLeaves.loci"

#intersectBed -a $MERGED/"CG-DMR-loci_collapsed.bed" -b $CG_LEAF/"CG-DMR-loci_collapsed.bed" -v |\
#awk '{OFS="\t"} {print $1, $2, $3}' > $OUTPATH/"CG-DMR.onlyFruits.loci"

## get the meth
#cat $CG_METH/"CG-DMR.merged.methylation.bed" |\
#sortBed -i - | intersectBed -a - -b $OUTPATH/"CG-DMR.common.loci" > $OUTPATH/"CG-DMR.common.loci.bed"
#cat $CG_METH/"CG-DMR.merged.methylation.bed" |\
#sortBed -i - | intersectBed -a - -b $OUTPATH/"CG-DMR.onlyLeaves.loci" > $OUTPATH/"CG-DMR.onlyLeaves.loci.bed"
#cat $CG_METH/"CG-DMR.merged.methylation.bed" | sortBed -i - |\
#intersectBed -a - -b $OUTPATH/"CG-DMR.onlyFruits.loci" > $OUTPATH/"CG-DMR.onlyFruits.loci.bed"

echo '======= C-DMR ======'
#intersectBed -a $C_LEAF/"C-DMR-loci_collapsed.bed" -b $C_FRUIT/"C-DMR-loci_collapsed.bed" |\
#awk '{OFS="\t"} {print $1, $2, $3}' > $OUTPATH/"C-DMR.common.loci"
#intersectBed -a $MERGED/"C-DMR-loci_collapsed.bed" -b $C_FRUIT/"C-DMR-loci_collapsed.bed" -v |\
#awk '{OFS="\t"} {print $1, $2, $3}' > $OUTPATH/"C-DMR.onlyLeaves.loci"
#intersectBed -a $MERGED/"C-DMR-loci_collapsed.bed" -b $C_LEAF/"C-DMR-loci_collapsed.bed" -v |\
#awk '{OFS="\t"} {print $1, $2, $3}' > $OUTPATH/"C-DMR.onlyFruits.loci"

## get the meth
cat $C_METH/"C-DMR.merged.methylation.bed" | sortBed -i - | intersectBed -a - -b $OUTPATH/"C-DMR.common.loci" > $OUTPATH/"C-DMR.common.loci.bed"
cat $C_METH/"C-DMR.merged.methylation.bed" | sortBed -i - | intersectBed -a - -b $OUTPATH/"C-DMR.onlyLeaves.loci" > $OUTPATH/"C-DMR.onlyLeaves.loci.bed"
cat $C_METH/"C-DMR.merged.methylation.bed" | sortBed -i - | intersectBed -a - -b $OUTPATH/"C-DMR.onlyFruits.loci" > $OUTPATH/"C-DMR.onlyFruits.loci.bed"
