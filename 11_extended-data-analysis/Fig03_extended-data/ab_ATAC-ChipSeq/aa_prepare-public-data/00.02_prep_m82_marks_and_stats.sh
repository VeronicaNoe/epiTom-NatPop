#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
### USER CONFIG
TAR="/mnt/disk2/vibanez/otherAnalysis/00_public-data/GSE245529_RAW.tar"
OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/aa_prepare-public-data"
M82_FASTA="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SLM82.fasta"
M82_GFF3_GZ="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SollycM82_genes_v1.1.1.gff3.gz"
# promoter upstream window (bp); change to 3000 if you want ±3kb sensitivity later
PROM_UP=3500
# Marks to extract + summarize
MARKS=(H3K9ac H3K27ac H3K9me2 H3K4me3 H3K36me3 Pol2 H3K27me3)
### CHECK TOOLS
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing tool: $1" >&2; exit 1; }; }
need tar
need awk
need sort
need bedtools
need samtools
need zcat
need grep

### DIRS
RAW_DIR="${OUTROOT}/atlas_peaks/raw"
MERGED_DIR="${OUTROOT}/atlas_peaks/merged"
MERGED_CHR_DIR="${OUTROOT}/atlas_peaks/merged_chr"
mkdir -p "$RAW_DIR" "$MERGED_DIR" "$MERGED_CHR_DIR" "$OUTROOT"

### 1) Build chrom sizes and canonical chr list (chr1-12)
samtools faidx "$M82_FASTA" >/dev/null 2>&1 || true
cut -f1,2 "${M82_FASTA}.fai" > "${OUTROOT}/M82.chrom.sizes"

awk '$1 ~ /^chr([1-9]|1[0-2])$/' "${OUTROOT}/M82.chrom.sizes" | cut -f1 > "${OUTROOT}/M82.canonical_chroms.txt"

GENOME_BP=$(awk 'FNR==NR{ok[$1]=1; next} ($1 in ok){s+=$2} END{print s}' \
  "${OUTROOT}/M82.canonical_chroms.txt" "${OUTROOT}/M82.chrom.sizes")

### 2) Build gene bodies BED and promoter BED (chr1-12 only)
# gene bodies
zcat "$M82_GFF3_GZ" \
| awk 'BEGIN{FS=OFS="\t"} $3=="gene"{ print $1, $4-1, $5 }' \
| sort -k1,1 -k2,2n \
| grep -F -f "${OUTROOT}/M82.canonical_chroms.txt" \
> "${OUTROOT}/M82.genes.chr.bed"

# promoters (strand-aware): [-PROM_UP..-1] from TSS
zcat "$M82_GFF3_GZ" \
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
}' "${OUTROOT}/M82.chrom.sizes" - \
| grep -F -f "${OUTROOT}/M82.canonical_chroms.txt" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> "${OUTROOT}/M82.promoters${PROM_UP}.chr.merged.bed"

### 3) Extract peak files for each mark from TAR (auto-discover)
# We accept either narrowPeak or broadPeak for a mark, because some are broad (e.g. H3K9me2, often H3K36me3).
echo "== Finding and extracting peak files from: $TAR"
for mark in "${MARKS[@]}"; do
  hits=$(tar -tf "$TAR" | grep -i -E "M82_${mark}.*peaks\.(narrowPeak|broadPeak)\.gz$" || true)
  if [[ -z "$hits" ]]; then
    echo "WARNING: no peak file found for $mark in TAR"
    continue
  fi

  echo "  $mark:"
  echo "$hits" | sed 's/^/    - /'

  # extract all replicates/files for that mark
  while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    tar -xf "$TAR" -C "$RAW_DIR" "$f"
  done <<< "$hits"
done

### 4) Union/merge peaks per mark -> union.bed and union.chr.bed
for mark in "${MARKS[@]}"; do
  files=( $(find "$RAW_DIR" -maxdepth 1 -type f | grep -i -E "M82_${mark}.*peaks\.(narrowPeak|broadPeak)\.gz$" || true) )
  if [[ ${#files[@]} -eq 0 ]]; then
    continue
  fi

  # union merge
  zcat -f "${files[@]}" \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "${MERGED_DIR}/${mark}.union.bed"

  # restrict to chr1-12
  grep -F -f "${OUTROOT}/M82.canonical_chroms.txt" "${MERGED_DIR}/${mark}.union.bed" \
  > "${MERGED_CHR_DIR}/${mark}.union.chr.bed"
done

### 5) Summary tables
STATS_TSV="${OUTROOT}/mark_stats.tsv"
CHR_TSV="${OUTROOT}/mark_chr_bp.tsv"

echo -e "mark\tn_peaks\tbp_covered\tpct_genome\tmedian_len\tpromoter_bp\tpromoter_pct_of_mark\tgene_bp\tgene_pct_of_mark" > "$STATS_TSV"
echo -e "mark\tchr\tbp_covered" > "$CHR_TSV"

for mark in "${MARKS[@]}"; do
  f="${MERGED_CHR_DIR}/${mark}.union.chr.bed"
  [[ -s "$f" ]] || continue

  n=$(wc -l < "$f")
  bp=$(awk '{s+=($3-$2)} END{print s+0}' "$f")
  pct=$(awk -v bp="$bp" -v G="$GENOME_BP" 'BEGIN{printf "%.3f",100*bp/G}')
  med=$(awk '{print $3-$2}' "$f" | sort -n | awk '{a[NR]=$1} END{if(NR%2) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}')

  prom_bp=$(bedtools intersect -a "$f" -b "${OUTROOT}/M82.promoters${PROM_UP}.chr.merged.bed" -wo \
            | awk '{s+=$NF} END{print s+0}')
  gene_bp=$(bedtools intersect -a "$f" -b "${OUTROOT}/M82.genes.chr.bed" -wo \
            | awk '{s+=$NF} END{print s+0}')

  prom_pct=$(awk -v pb="$prom_bp" -v tb="$bp" 'BEGIN{printf "%.2f", (tb>0?100*pb/tb:0)}')
  gene_pct=$(awk -v gb="$gene_bp" -v tb="$bp" 'BEGIN{printf "%.2f", (tb>0?100*gb/tb:0)}')

  echo -e "${mark}\t${n}\t${bp}\t${pct}\t${med}\t${prom_bp}\t${prom_pct}\t${gene_bp}\t${gene_pct}" >> "$STATS_TSV"

  # chromosome bp distribution
  awk -v m="$mark" 'BEGIN{OFS="\t"}{bp[$1]+=($3-$2)} END{for(c in bp) print m,c,bp[c]}' "$f" \
  | sort -k2,2 -k3,3nr >> "$CHR_TSV"
done

echo "Wrote:"
echo "  $STATS_TSV"
echo "  $CHR_TSV"
echo "Peaks:"
echo "  raw:        $RAW_DIR"
echo "  merged:     $MERGED_DIR"
echo "  merged_chr: $MERGED_CHR_DIR"
