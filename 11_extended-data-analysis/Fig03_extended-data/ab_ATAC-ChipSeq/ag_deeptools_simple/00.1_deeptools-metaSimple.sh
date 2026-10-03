#!/usr/bin/env bash
DMRSET="$( cat $1 | cut -d'_' -f1)"
REG="$( cat $1 | cut -d'_' -f2)"

OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq"
BW="${OUTROOT}/aa_prepare-public-data/atlas_bw"
OUTDT="${OUTROOT}/ag_deeptools_simple"

NCPU="${NCPU:-8}"

# -------------------------
# Resolve input BEDs
# -------------------------
# Low/high fUM bins (exported by your R step)
# Prefer the "ae_deeptools_out/deeptools_beds" layout, fallback to OUTROOT/deeptools_beds.
BEDROOT1="${OUTROOT}/ae_deeptools_out/deeptools_beds"
BEDROOT2="${OUTROOT}/deeptools_beds"
if [[ -d "$BEDROOT1" ]]; then
  BEDROOT="$BEDROOT1"
elif [[ -d "$BEDROOT2" ]]; then
  BEDROOT="$BEDROOT2"
else
  echo "ERROR: cannot find deeptools_beds in $BEDROOT1 or $BEDROOT2" >&2
  exit 2
fi

LOW="${BEDROOT}/${DMRSET}/${DMRSET}.${REG}.fUM_0p0-0p1.bed"
HIGH="${BEDROOT}/${DMRSET}/${DMRSET}.${REG}.fUM_0p9-1p0.bed"
if [[ ! -s "$LOW" || ! -s "$HIGH" ]]; then
  echo "WARNING: missing LOW/HIGH beds for ${DMRSET} ${REG} (skip)" >&2
  exit 0
fi

# Region universe + genome sizes (M82 chr1-12)
REGROOT="${OUTROOT}/ac_annotate-DMRs-M82/regions_M82"
GENOME_CHR="${REGROOT}/M82.chrom.sizes.chr"

UNIV_PROM="${REGROOT}/promoter3500.merged.bed3"
UNIV_GENE="${REGROOT}/genes.merged.bed3"
UNIV_GENE_TE="${REGROOT}/genes_TE.bed3"
UNIV_TE="${REGROOT}/TE.merged.bed3"
UNIV_INTER="${REGROOT}/intergenic.bed3"

get_universe () {
  case "$1" in
    promoter)   echo "$UNIV_PROM" ;;
    gene)       echo "$UNIV_GENE" ;;
    gene-TE|gene_TE) echo "$UNIV_GENE_TE" ;;
    TE|tes)     echo "$UNIV_TE" ;;
    intergenic) echo "$UNIV_INTER" ;;
    *)          echo "" ;;
  esac
}
UNIV="$(get_universe "$REG")"
if [[ -z "$UNIV" || ! -s "$UNIV" ]]; then
  echo "ERROR: universe BED missing/unknown REG=$REG: $UNIV" >&2
  exit 2
fi
[[ -s "$GENOME_CHR" ]] || { echo "ERROR: missing genome file $GENOME_CHR" >&2; exit 2; }

# Region-wide DMRs (all DMRs in that region) -> used to generate ONE shuffle per region
REG_DMR_BED4="${OUTROOT}/ac_annotate-DMRs-M82/dmrs_by_region/${DMRSET}.${REG}.bed4"
if [[ ! -s "$REG_DMR_BED4" ]]; then
  echo "WARNING: missing region DMR bed4: $REG_DMR_BED4 (skip)" >&2
  exit 0
fi

# -------------------------
# Signals: ChIP + Pol2 + ATAC
# -------------------------
SIGS=()
SLABELS=()

check_bw () {
  local f="$1"
  python - <<PY 2>/dev/null
import pyBigWig
bw = pyBigWig.open("${f}")
bw.close()
PY
}

add_sig () {
  local f="$1"
  local lab="$2"
  if [[ -n "${f}" && -s "${f}" ]]; then
    if check_bw "$f"; then
      SIGS+=("$f")
      SLABELS+=("$lab")
    else
      echo "WARNING: bigWig unreadable (skip): $f" >&2
    fi
  else
    echo "WARNING: missing bigWig for ${lab} (skip)" >&2
  fi
}

H3K4ME3=$(ls -1 "${BW}"/*H3K4me3*.bigwig 2>/dev/null | head -n1 || true)
H3K36ME3=$(ls -1 "${BW}"/*H3K36me3*.bigwig 2>/dev/null | head -n1 || true)
H3K27AC=$(ls -1 "${BW}"/*H3K27ac*.bigwig 2>/dev/null | head -n1 || true)
H3K9AC=$(ls -1 "${BW}"/*H3K9ac*.bigwig 2>/dev/null | head -n1 || true)
H3K9ME2=$(ls -1 "${BW}"/*H3K9me2*.bigwig 2>/dev/null | head -n1 || true)
POL2=$(ls -1 "${BW}"/*Pol2*.bigwig 2>/dev/null | head -n1 || true)

ATAC1=$(ls -1 "${BW}"/*ATAC*M82*rep1*.bigwig 2>/dev/null | head -n1 || true)
ATAC2=$(ls -1 "${BW}"/*ATAC*M82*rep2*.bigwig 2>/dev/null | head -n1 || true)
if [[ -z "$ATAC1" ]]; then ATAC1=$(ls -1 "${BW}"/*ATAC*M82*.bigwig 2>/dev/null | head -n1 || true); fi
if [[ -z "$ATAC2" ]]; then ATAC2=$(ls -1 "${BW}"/*ATAC*M82*.bigwig 2>/dev/null | sed -n '2p' || true); fi

add_sig "$H3K4ME3" "H3K4me3"
add_sig "$H3K36ME3" "H3K36me3"
add_sig "$H3K27AC"  "H3K27ac"
add_sig "$H3K9AC"   "H3K9ac"
add_sig "$POL2"     "Pol2"
add_sig "$ATAC1"    "ATAC_rep1"
add_sig "$ATAC2"    "ATAC_rep2"
add_sig "$H3K9ME2"  "H3K9me2"

if [[ ${#SIGS[@]} -eq 0 ]]; then
  echo "ERROR: no readable bigWigs found in ${BW}" >&2
  exit 1
fi

# -------------------------
# Build ONE region-wide shuffle (lengths match the region DMRs)
# -------------------------
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

REG_BED3="${TMPDIR}/${DMRSET}.${REG}.bed3"
cut -f1-3 "$REG_DMR_BED4" > "$REG_BED3"

SHUF="${OUTDT}/${DMRSET}.${REG}.regionShuffle.bed3"
bedtools shuffle -i "$REG_BED3" -g "$GENOME_CHR" -incl "$UNIV" -chrom -seed 123 \
| bedtools sort -g "$GENOME_CHR" -i - \
> "$SHUF"

NLOW=$(wc -l < "$LOW" | awk '{print $1}')
NHIGH=$(wc -l < "$HIGH" | awk '{print $1}')
NSHUF=$(wc -l < "$SHUF" | awk '{print $1}')

# -------------------------
# deepTools: 3 groups in one PDF (low, high, shuffle-region)
# referencePoint center = midpoint of each DMR interval
# -------------------------
MATRIX="${OUTDT}/${DMRSET}.${REG}.low_high_vs_regionShuffle.matrix.gz"
PDF="${OUTDT}/${DMRSET}.${REG}.low_high_vs_regionShuffle.profile.pdf"

computeMatrix reference-point \
  --referencePoint center \
  -b 2000 -a 2000 \
  -R "$LOW" "$HIGH" "$SHUF" \
  -S "${SIGS[@]}" \
  --skipZeros \
  --numberOfProcessors "$NCPU" \
  -o "$MATRIX"

plotProfile -m "$MATRIX" -out "$PDF" --perGroup \
  --regionsLabel "low fUM (n=${NLOW})" "high fUM (n=${NHIGH})" "region shuffle (n=${NSHUF})" \
  --samplesLabel "${SLABELS[@]}"

echo "Wrote: $PDF"
