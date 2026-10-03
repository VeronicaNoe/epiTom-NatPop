#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

DMRSET="$(cat "$1")"

OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/ac_annotate-DMRs-M82"
INDIR="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/aa_prepare-public-data"
LIFTROOT="/mnt/disk2/vibanez/otherAnalysis/04_liftoff-DMRs-SL2.5-M82"
GENEPREP="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation"

DMR_BED4_GZ="${LIFTROOT}/${DMRSET}.M82.hiconf.mapq30.chr.bed4.gz"

PEAKDIR="${INDIR}/atlas_peaks/merged_chr"
MARKS=(H3K9ac H3K27ac H3K9me2 H3K4me3 H3K36me3 Pol2 H3K27me3)

# Your promoter/gene beds (already built)
GENES_BED="${INDIR}/M82.genes.chr.bed"
PROM_BED="${INDIR}/M82.promoters3500.chr.merged.bed"

# Genome order file (chr1–chr12)
CHROMSIZES_CHR="${OUTROOT}/regions_M82/M82.chrom.sizes.chr"
CANON_CHR="${OUTROOT}/regions_M82/M82.canonical_chroms.txt"

# TE annotation input (EDTA)
TE_GFF3="${GENEPREP}/SLM82.fa.mod.EDTA.TEanno.gff3"
TE_BED_MERGED="${OUTROOT}/regions_M82/M82.TE.chr.merged.bed3"

OUTDIR="${OUTROOT}/dmr_mark_annotations/${DMRSET}"
mkdir -p "$OUTDIR"

DMR_BED4="${OUTDIR}/${DMRSET}.bed4"
BASE_TSV="${OUTDIR}/${DMRSET}.base.tsv"
OUT_TSV_GZ="${OUTDIR}/${DMRSET}.dmr_chromatin_annot.tsv.gz"

# -----------------------------
# helpers
# -----------------------------
cov_cols () {
  local A="$1"
  local B="$2"
  bedtools coverage -a "$A" -b "$B" \
  | awk 'BEGIN{FS=OFS="\t"}{ bp=$(NF-2)+0; frac=$NF+0; any=(bp>0?1:0); print any,bp,frac }'
}

# -----------------------------
# build TE merged bed if missing
# -----------------------------
if [[ ! -s "$TE_BED_MERGED" ]]; then
  echo "Building TE BED (chr1-12) from: $TE_GFF3"
  mkdir -p "$(dirname "$TE_BED_MERGED")"

  # If canonical chr list is missing, derive it from chromsizes
  if [[ ! -s "$CANON_CHR" ]]; then
    cut -f1 "$CHROMSIZES_CHR" > "$CANON_CHR"
  fi

  # Convert GFF3 to BED3 and merge, respecting genome order
  # Note: EDTA gff3 may contain various feature types; we take all non-comment lines.
  awk 'BEGIN{FS=OFS="\t"} $0!~/^#/{print $1,$4-1,$5}' "$TE_GFF3" \
  | grep -F -f "$CANON_CHR" \
  | bedtools sort -g "$CHROMSIZES_CHR" -i - \
  | bedtools merge -i - \
  > "$TE_BED_MERGED"

  echo "Wrote: $TE_BED_MERGED"
fi

# -----------------------------
# prep DMR base
# -----------------------------
if [[ ! -s "$DMR_BED4_GZ" ]]; then
  echo "ERROR: missing $DMR_BED4_GZ" >&2
  exit 1
fi

zcat "$DMR_BED4_GZ" > "$DMR_BED4"
cut -f1-4 "$DMR_BED4" > "$BASE_TSV"

# -----------------------------
# compute columns
# -----------------------------
COLFILES=()

# marks
for m in "${MARKS[@]}"; do
  f="${PEAKDIR}/${m}.union.chr.bed"
  if [[ ! -s "$f" ]]; then
    echo "WARNING: missing $f (skip $m)" >&2
    continue
  fi
  out="${OUTDIR}/${m}.cols.tsv"
  cov_cols "$DMR_BED4" "$f" > "$out"
  COLFILES+=("$out")
done

# promoters
if [[ -s "$PROM_BED" ]]; then
  out="${OUTDIR}/PROMOTER.cols.tsv"
  cov_cols "$DMR_BED4" "$PROM_BED" > "$out"
  COLFILES+=("$out")
else
  echo "WARNING: missing PROM_BED=$PROM_BED (skip promoters)" >&2
fi

# genes
if [[ -s "$GENES_BED" ]]; then
  out="${OUTDIR}/GENE.cols.tsv"
  cov_cols "$DMR_BED4" "$GENES_BED" > "$out"
  COLFILES+=("$out")
else
  echo "WARNING: missing GENES_BED=$GENES_BED (skip genes)" >&2
fi

# TEs
if [[ -s "$TE_BED_MERGED" ]]; then
  out="${OUTDIR}/TE.cols.tsv"
  cov_cols "$DMR_BED4" "$TE_BED_MERGED" > "$out"
  COLFILES+=("$out")
else
  echo "WARNING: missing TE_BED_MERGED=$TE_BED_MERGED (skip TE)" >&2
fi

# -----------------------------
# header
# -----------------------------
{
  printf "chr\tstart\tend\tid"
  for m in "${MARKS[@]}"; do
    f="${PEAKDIR}/${m}.union.chr.bed"
    [[ -s "$f" ]] || continue
    printf "\t%s_any\t%s_bp\t%s_frac" "$m" "$m" "$m"
  done
  [[ -s "$PROM_BED" ]] && printf "\tPROMOTER_any\tPROMOTER_bp\tPROMOTER_frac"
  [[ -s "$GENES_BED" ]] && printf "\tGENE_any\tGENE_bp\tGENE_frac"
  [[ -s "$TE_BED_MERGED" ]] && printf "\tTE_any\tTE_bp\tTE_frac"
  printf "\n"
} > "${OUTDIR}/header.tsv"

# -----------------------------
# write output
# -----------------------------
paste "$BASE_TSV" "${COLFILES[@]}" \
| cat "${OUTDIR}/header.tsv" - \
| gzip -c > "$OUT_TSV_GZ"

echo "Wrote: $OUT_TSV_GZ"
