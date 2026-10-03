DIR="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results"
ANNO_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"

sort -k1,1 -k2,2n "${ANNO_DIR}/allGeneNames-Function.SL5.sorted.bed" > "${ANNO_DIR}/allGeneNames-Function.sorted.bed"

awk -F',' 'BEGIN{OFS="\t"} NR>1 {
  chr=sprintf("%02d",$2)
  pos=$3
  print chr, pos-1, pos, $0
}' "${DIR}/02.2_CG-DMR_top_1_percent.csv" \
| sort -k1,1 -k2,2n \
| bedtools closest -a - -b "${ANNO_DIR}/allGeneNames-Function.sorted.bed" -d \
| awk 'BEGIN{OFS="\t"} $NF <= 10000' \
> "${DIR}/02.3_CG-DMR.closest-gene_top-1-SVs.tsv"

awk -F',' 'BEGIN{OFS="\t"} NR>1 {
  chr=sprintf("%02d",$2)
  pos=$3
  print chr, pos-1, pos, $0
}' "${DIR}/02.2_C-DMR_top_1_percent.csv" \
| sort -k1,1 -k2,2n \
| bedtools closest -a - -b "${ANNO_DIR}/allGeneNames-Function.sorted.bed" -d \
| awk 'BEGIN{OFS="\t"} $NF <= 10000' \
> "${DIR}/02.3_C-DMR.closest-gene_top-1-SVs.tsv"

