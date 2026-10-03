#!/usr/bin/env bash
REG="$( cat $1 | cut -d'_' -f1)"
OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq"
BW="${OUTROOT}/aa_prepare-public-data/atlas_bw"
OUTDT="/mnt/disk2/vibanez/otherAnalysis/review_chr-marks"

NCPU="${NCPU:-8}"
UP="${UP:-3000}"
DOWN="${DOWN:-3000}"
SEED="${SEED:-123}"

mkdir -p "$OUTDT"

# -------------------------
# Input BED roots
# -------------------------
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

# -------------------------
# Region universes + genome
# -------------------------
REGROOT="${OUTROOT}/ac_annotate-DMRs-M82/regions_M82"
GENOME_CHR="${REGROOT}/M82.chrom.sizes.chr"

UNIV_PROM="${REGROOT}/promoter3500.merged.bed3"
UNIV_GENE="${REGROOT}/genes.merged.bed3"
UNIV_GENE_TE="${REGROOT}/genes_TE.bed3"
UNIV_TE="${REGROOT}/TE.merged.bed3"
UNIV_INTER="${REGROOT}/intergenic.bed3"

get_universe () {
  case "$1" in
    promoter)           echo "$UNIV_PROM" ;;
    gene)               echo "$UNIV_GENE" ;;
    gene-TE|gene_TE)    echo "$UNIV_GENE_TE" ;;
    TE|tes)             echo "$UNIV_TE" ;;
    intergenic)         echo "$UNIV_INTER" ;;
    *)
      echo ""
      ;;
  esac
}

UNIV="$(get_universe "$REG")"
if [[ -z "$UNIV" || ! -s "$UNIV" ]]; then
  echo "ERROR: universe BED missing/unknown REG=$REG : $UNIV" >&2
  exit 2
fi
[[ -s "$GENOME_CHR" ]] || { echo "ERROR: missing genome file $GENOME_CHR" >&2; exit 2; }

# -------------------------
# DMR classes
# -------------------------
DMRSETS=("CG-DMR" "C-DMR")

# -------------------------
# Resolve BEDs
# -------------------------
declare -A LOWBED HIGHBED SHUFBED NLOW NHIGH NSHUF

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

for DMRSET in "${DMRSETS[@]}"; do
  LOW="${BEDROOT}/${DMRSET}/${DMRSET}.${REG}.fUM_0p0-0p1.bed"
  HIGH="${BEDROOT}/${DMRSET}/${DMRSET}.${REG}.fUM_0p9-1p0.bed"

  if [[ ! -s "$LOW" ]]; then
    echo "ERROR: missing LOW bed: $LOW" >&2
    exit 2
  fi
  if [[ ! -s "$HIGH" ]]; then
    echo "ERROR: missing HIGH bed: $HIGH" >&2
    exit 2
  fi

  REG_DMR_BED4="${OUTROOT}/ac_annotate-DMRs-M82/dmrs_by_region/${DMRSET}.${REG}.bed4"
  if [[ ! -s "$REG_DMR_BED4" ]]; then
    echo "ERROR: missing region DMR bed4: $REG_DMR_BED4" >&2
    exit 2
  fi

  REG_BED3="${TMPDIR}/${DMRSET}.${REG}.bed3"
  cut -f1-3 "$REG_DMR_BED4" > "$REG_BED3"

  SHUF="${OUTDT}/${DMRSET}.${REG}.regionShuffle.bed3"
  bedtools shuffle \
    -i "$REG_BED3" \
    -g "$GENOME_CHR" \
    -incl "$UNIV" \
    -chrom \
    -seed "$SEED" \
  | bedtools sort -g "$GENOME_CHR" -i - \
  > "$SHUF"

  LOWBED["$DMRSET"]="$LOW"
  HIGHBED["$DMRSET"]="$HIGH"
  SHUFBED["$DMRSET"]="$SHUF"

  NLOW["$DMRSET"]="$(wc -l < "$LOW" | awk '{print $1}')"
  NHIGH["$DMRSET"]="$(wc -l < "$HIGH" | awk '{print $1}')"
  NSHUF["$DMRSET"]="$(wc -l < "$SHUF" | awk '{print $1}')"
done

# -------------------------
# Signals
# Recommended order:
# ATAC -> Pol2 -> active marks -> repressive mark
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
  if [[ -n "$f" && -s "$f" ]]; then
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

ATAC1=$(ls -1 "${BW}"/*ATAC*M82*rep1*.bigwig 2>/dev/null | head -n1 || true)
ATAC2=$(ls -1 "${BW}"/*ATAC*M82*rep2*.bigwig 2>/dev/null | head -n1 || true)
if [[ -z "$ATAC1" ]]; then ATAC1=$(ls -1 "${BW}"/*ATAC*M82*.bigwig 2>/dev/null | head -n1 || true); fi
if [[ -z "$ATAC2" ]]; then ATAC2=$(ls -1 "${BW}"/*ATAC*M82*.bigwig 2>/dev/null | sed -n '2p' || true); fi

POL2=$(ls -1 "${BW}"/*Pol2*.bigwig 2>/dev/null | head -n1 || true)
H3K27AC=$(ls -1 "${BW}"/*H3K27ac*.bigwig 2>/dev/null | head -n1 || true)
H3K9AC=$(ls -1 "${BW}"/*H3K9ac*.bigwig 2>/dev/null | head -n1 || true)
H3K4ME3=$(ls -1 "${BW}"/*H3K4me3*.bigwig 2>/dev/null | head -n1 || true)
H3K36ME3=$(ls -1 "${BW}"/*H3K36me3*.bigwig 2>/dev/null | head -n1 || true)
H3K9ME2=$(ls -1 "${BW}"/*H3K9me2*.bigwig 2>/dev/null | head -n1 || true)

add_sig "$ATAC1"    "ATAC_rep1"
add_sig "$ATAC2"    "ATAC_rep2"
add_sig "$POL2"     "Pol2"
add_sig "$H3K27AC"  "H3K27ac"
add_sig "$H3K9AC"   "H3K9ac"
add_sig "$H3K4ME3"  "H3K4me3"
add_sig "$H3K36ME3" "H3K36me3"
add_sig "$H3K9ME2"  "H3K9me2"

if [[ ${#SIGS[@]} -eq 0 ]]; then
  echo "ERROR: no readable bigWigs found in ${BW}" >&2
  exit 1
fi

# -------------------------
# One matrix for this annotation:
# CG low/high/shuf + C low/high/shuf
# -------------------------
MATRIX="${OUTDT}/${REG}.CG_C.low_high_shuffle.matrix.gz"
PDF="${OUTDT}/${REG}.CG_C.low_high_shuffle.profile.pdf"

computeMatrix reference-point \
  --referencePoint center \
  -b "$UP" -a "$DOWN" \
  -R \
    "${LOWBED[CG-DMR]}" \
    "${HIGHBED[CG-DMR]}" \
    "${SHUFBED[CG-DMR]}" \
    "${LOWBED[C-DMR]}" \
    "${HIGHBED[C-DMR]}" \
    "${SHUFBED[C-DMR]}" \
  -S "${SIGS[@]}" \
  --skipZeros \
  --numberOfProcessors "$NCPU" \
  -o "$MATRIX"

plotProfile \
  -m "$MATRIX" \
  -out "$PDF" \
  --perGroup \
  --regionsLabel \
    "CG-DMR low fUM (n=${NLOW[CG-DMR]})" \
    "CG-DMR high fUM (n=${NHIGH[CG-DMR]})" \
    "CG-DMR region shuffle (n=${NSHUF[CG-DMR]})" \
    "C-DMR low fUM (n=${NLOW[C-DMR]})" \
    "C-DMR high fUM (n=${NHIGH[C-DMR]})" \
    "C-DMR region shuffle (n=${NSHUF[C-DMR]})" \
  --samplesLabel "${SLABELS[@]}"

echo "Wrote:"
echo "  $MATRIX"
echo "  $PDF"
