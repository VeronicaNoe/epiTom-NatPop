#!/bin/bash
KO="$( cat $1 )"
echo $KO
KO_DIR="/mnt/disk2/vibanez/natural-experimental_DMR/ab_KO"
#"/mnt/disk2/vibanez/09_KO-processing/09.0_biseq/ae_DMR/ba_DMR-calling"
DMR_DIR="/mnt/disk2/vibanez/natural-experimental_DMR/aa_DMR"
OUT_DIR="/mnt/disk2/vibanez/natural-experimental_DMR/ac_output"

#C
cat $KO_DIR/$KO"_CG-DMR.methylation.bed" | awk -v KO=$KO '{OFS="\t"}{print $1,$2,$3, $4=KO}' | intersectBed -a - -b $DMR_DIR/"C-DMR.bed" -nonamecheck > $OUT_DIR/$KO"_CG_C-DMR.tmp"
cat $KO_DIR/$KO"_CHG-DMR.methylation.bed" | awk -v KO=$KO '{OFS="\t"}{print $1,$2,$3, $4=KO}'| intersectBed -a - -b $DMR_DIR/"C-DMR.bed" -nonamecheck > $OUT_DIR/$KO"_CHG_C-DMR.tmp"
cat $KO_DIR/$KO"_CHH-DMR.methylation.bed" | awk -v KO=$KO '{OFS="\t"}{print $1,$2,$3, $4=KO}'| intersectBed -a - -b $DMR_DIR/"C-DMR.bed" -nonamecheck > $OUT_DIR/$KO"_CHH_C-DMR.tmp"

cat $OUT_DIR/$KO*"_C-DMR.tmp" | sort -k1,1 -k2,2n | mergeBed -i - -c 4 -o distinct > $OUT_DIR/$KO".C-DMR.toMerge"
cat $KO_DIR/$KO*".methylation.bed" | sort -k1,1 -k2,2n | mergeBed -i - | intersectBed -a - -b $OUT_DIR/$KO".C-DMR.toMerge" -v > $OUT_DIR/$KO".C-DMR.only"

cat $DMR_DIR/"C-DMR.bed" | awk '{OFS="\t"}{print $1,$2,$3, $4=1}' | intersectBed -a - -b $KO_DIR/$KO"_CG-DMR.methylation.bed" -nonamecheck -v > $OUT_DIR/$KO"_CG_nat_C-DMR.tmp"
cat $DMR_DIR/"C-DMR.bed" | awk '{OFS="\t"}{print $1,$2,$3, $4=1}' | intersectBed -a - -b $KO_DIR/$KO"_CHG-DMR.methylation.bed" -nonamecheck -v > $OUT_DIR/$KO"_CHG_nat_C-DMR.tmp"
cat $DMR_DIR/"C-DMR.bed" | awk '{OFS="\t"}{print $1,$2,$3, $4=1}' | intersectBed -a - -b $KO_DIR/$KO"_CHH-DMR.methylation.bed" -nonamecheck -v > $OUT_DIR/$KO"_CHH_nat_C-DMR.tmp"

cat $OUT_DIR/$KO*"_nat_C-DMR.tmp" | sort -k1,1 -k2,2n | mergeBed -i - -c 4 -o distinct > $OUT_DIR/$KO"_nat.C-DMR.toMerge"

cat $OUT_DIR/$KO*"toMerge" |  sort -k1,1 -k2,2n | mergeBed -i - -c 4 -o distinct > $OUT_DIR/$KO"_C-DMR.bed"
rm $OUT_DIR/$KO*".tmp"
rm $OUT_DIR/$KO*".toMerge"


#CG
cat $KO_DIR/$KO"_CG-DMR.methylation.bed" | awk -v KO=$KO '{OFS="\t"}{print $1,$2,$3, $4=KO}' | intersectBed -a - -b $DMR_DIR/"CG-DMR.bed" -nonamecheck > $OUT_DIR/$KO"_CG_CG-DMR.tmp"
cat $KO_DIR/$KO"_CHG-DMR.methylation.bed" | awk -v KO=$KO '{OFS="\t"}{print $1,$2,$3, $4=KO}' | intersectBed -a - -b $DMR_DIR/"CG-DMR.bed" -nonamecheck > $OUT_DIR/$KO"_CHG_CG-DMR.tmp"
cat $KO_DIR/$KO"_CHH-DMR.methylation.bed" | awk -v KO=$KO '{OFS="\t"}{print $1,$2,$3, $4=KO}' | intersectBed -a - -b $DMR_DIR/"CG-DMR.bed" -nonamecheck > $OUT_DIR/$KO"_CHH_CG-DMR.tmp"

cat $OUT_DIR/$KO*"_CG-DMR.tmp" | sort -k1,1 -k2,2n | mergeBed -i - -c 4 -o distinct > $OUT_DIR/$KO".CG-DMR.toMerge"
cat $KO_DIR/$KO*".methylation.bed" | sort -k1,1 -k2,2n | mergeBed -i - | intersectBed -a - -b $OUT_DIR/$KO".CG-DMR.toMerge" -v > $OUT_DIR/$KO".CG-DMR.only"

cat $DMR_DIR/"CG-DMR.bed" | awk '{OFS="\t"}{print $1,$2,$3, $4=1}' | intersectBed -a - -b $KO_DIR/$KO"_CG-DMR.methylation.bed" -nonamecheck -v > $OUT_DIR/$KO"_CG_nat_CG-DMR.tmp"
cat $DMR_DIR/"CG-DMR.bed" | awk '{OFS="\t"}{print $1,$2,$3, $4=1}' | intersectBed -a - -b $KO_DIR/$KO"_CHG-DMR.methylation.bed" -nonamecheck -v > $OUT_DIR/$KO"_CHG_nat_CG-DMR.tmp"
cat $DMR_DIR/"CG-DMR.bed" | awk '{OFS="\t"}{print $1,$2,$3, $4=1}' | intersectBed -a - -b $KO_DIR/$KO"_CHH-DMR.methylation.bed" -nonamecheck -v > $OUT_DIR/$KO"_CHH_nat_CG-DMR.tmp"

cat $OUT_DIR/$KO*"_nat_CG-DMR.tmp" | sort -k1,1 -k2,2n | mergeBed -i - -c 4 -o distinct > $OUT_DIR/$KO"_nat.CG-DMR.toMerge"

cat $OUT_DIR/$KO*"toMerge" |  sort -k1,1 -k2,2n | mergeBed -i - -c 4 -o distinct > $OUT_DIR/$KO"_CG-DMR.bed"
rm $OUT_DIR/$KO*".tmp"
rm $OUT_DIR/$KO*".toMerge"

cat $OUT_DIR/$KO*".only" | sort -k1,1 -k2,2n | mergeBed -i - > $OUT_DIR/$KO"_only.bed"
rm $OUT_DIR/$KO*".only"
