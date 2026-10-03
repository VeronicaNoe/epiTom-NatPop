#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/ac_annotate-DMRs-M82"
GENEPREP="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation"
M82_FAI="${GENEPREP}/SLM82.fasta.fai"
GENES_GFF3_GZ="${GENEPREP}/SollycM82_genes_v1.1.1.gff3.gz"
TE_GFF3="${GENEPREP}/SLM82.fa.mod.EDTA.TEanno.gff3"
G="${OUTROOT}/regions_M82/M82.chrom.sizes.chr"

PROM_UP=3500

mkdir -p "${OUTROOT}/regions_M82"

# chrom sizes + canonical chr1-12
cut -f1,2 "$M82_FAI" > "${OUTROOT}/regions_M82/M82.chrom.sizes"
awk '$1 ~ /^chr([1-9]|1[0-2])$/' "${OUTROOT}/regions_M82/M82.chrom.sizes" > "${OUTROOT}/regions_M82/M82.chrom.sizes.chr"
cut -f1 "${OUTROOT}/regions_M82/M82.chrom.sizes.chr" > "${OUTROOT}/regions_M82/M82.canonical_chroms.txt"

# genes BED (chr only)
zcat "$GENES_GFF3_GZ" \
| awk 'BEGIN{FS=OFS="\t"} $3=="gene"{print $1,$4-1,$5,$7}' \
| grep -F -f "${OUTROOT}/regions_M82/M82.canonical_chroms.txt" \
| sort -k1,1 -k2,2n \
> "${OUTROOT}/regions_M82/genes.bed4"

cut -f1-3 "${OUTROOT}/regions_M82/genes.bed4" \
| bedtools sort -i - \
| bedtools merge -i - \
> "${OUTROOT}/regions_M82/genes.merged.bed3"

# promoters (strand-aware), output BED3 merged
# Uses chr lengths to clamp ends
awk 'BEGIN{FS=OFS="\t"} {L[$1]=$2} END{}' "${OUTROOT}/regions_M82/M82.chrom.sizes.chr" >/dev/null

zcat "$GENES_GFF3_GZ" \
| awk -v UP="$PROM_UP" 'BEGIN{FS=OFS="\t"}
FNR==NR{L[$1]=$2; next}
$3=="gene"{
  chr=$1; start=$4; end=$5; strand=$7;
  if(!(chr in L)) next;

  if(strand=="+"){
    p1 = start-UP; if(p1<1) p1=1;
    p2 = start-1;  if(p2<1) next;
    print chr, p1-1, p2;
  } else if(strand=="-"){
    p1 = end+1;
    p2 = end+UP; if(p2>L[chr]) p2=L[chr];
    if(p1>p2) next;
    print chr, p1-1, p2;
  }
}' "${OUTROOT}/regions_M82/M82.chrom.sizes.chr" - \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> "${OUTROOT}/regions_M82/promoter${PROM_UP}.merged.bed3"

# TE BED3 merged from EDTA GFF3 (chr only)
# (EDTA TEanno.gff3 is already TE features; we take all non-comment lines)
awk 'BEGIN{FS=OFS="\t"} $0!~/^#/{print $1,$4-1,$5}' "$TE_GFF3" \
| grep -F -f "${OUTROOT}/regions_M82/M82.canonical_chroms.txt" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> "${OUTROOT}/regions_M82/TE.merged.bed3"

# genes_TE vs genes_noTE (gene intervals)
bedtools intersect -u -a "${OUTROOT}/regions_M82/genes.merged.bed3" -b "${OUTROOT}/regions_M82/TE.merged.bed3" \
> "${OUTROOT}/regions_M82/genes_TE.bed3"

bedtools intersect -v -a "${OUTROOT}/regions_M82/genes.merged.bed3" -b "${OUTROOT}/regions_M82/TE.merged.bed3" \
> "${OUTROOT}/regions_M82/genes_noTE.bed3"

# intergenic universe = complement of (genes ∪ promoters ∪ TE)
cat "${OUTROOT}/regions_M82/genes.merged.bed3" \
    "${OUTROOT}/regions_M82/promoter${PROM_UP}.merged.bed3" \
    "${OUTROOT}/regions_M82/TE.merged.bed3" \
| bedtools sort -g "$G" -i - \
| bedtools merge -i - \
> "${OUTROOT}/regions_M82/union_genes_prom_TE.bed3"

bedtools complement \
  -i "${OUTROOT}/regions_M82/union_genes_prom_TE.bed3" \
  -g "${OUTROOT}/regions_M82/M82.chrom.sizes.chr" \
> "${OUTROOT}/regions_M82/intergenic.bed3"

echo "Wrote universes in: ${OUTROOT}/regions_M82"
ls -lh "${OUTROOT}/regions_M82" | sed -n '1,20p'

