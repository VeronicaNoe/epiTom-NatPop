#!/usr/bin/env bash
INDIR="${1:-/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/ag_deeptools_simple}"
OUTDIR="${2:-/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/ag_deeptools_simple}"
TAG="${3:-low_high_vs_regionShuffle}"
NCPU="${NCPU:-8}"

DMRSETS=(C-DMR CG-DMR)
REGIONS=(promoter gene gene-TE TE intergenic)

need_cmd () { command -v "$1" >/dev/null 2>&1; }
have_pdfjam=0
have_pdfunite=0
need_cmd pdfjam && have_pdfjam=1
need_cmd pdfunite && have_pdfunite=1

mkdir -p "$OUTDIR"

# resolve matrix path (gene-TE vs gene_TE variants)
matrix_path () {
  local d="$1" r="$2"
  local p1="${INDIR}/${d}.${r}.${TAG}.matrix.gz"
  local p2="${INDIR}/${d}.${r//-/_}.${TAG}.matrix.gz"
  if [[ -s "$p1" ]]; then
    echo "$p1"
  elif [[ -s "$p2" ]]; then
    echo "$p2"
  else
    echo ""
  fi
}

# compute max of the *mean profile* across groups/samples/bins for a matrix
max_mean_profile_from_matrix () {
  local mfile="$1"
  python - "$mfile" <<'PY'
import sys, numpy as np
from deeptools import heatmapper

mfile = sys.argv[1]
hm = heatmapper.heatmapper()
hm.read_matrix_file(mfile)
M = hm.matrix.matrix  # shape: (regions, bins*samples)

# boundaries are usually lists of ints of length (nGroups+1)/(nSamples+1)
gb = getattr(hm.matrix, "group_boundaries", None)
sb = getattr(hm.matrix, "sample_boundaries", None)
if gb is None or sb is None:
    raise SystemExit("Could not find group_boundaries/sample_boundaries in matrix object")

# normalize boundaries into slices
if len(gb) and isinstance(gb[0], (list, tuple)) and len(gb[0]) == 2:
    group_slices = gb
else:
    group_slices = list(zip(gb[:-1], gb[1:]))

if len(sb) and isinstance(sb[0], (list, tuple)) and len(sb[0]) == 2:
    sample_slices = sb
else:
    sample_slices = list(zip(sb[:-1], sb[1:]))

maxv = 0.0
for g0, g1 in group_slices:
    rows = M[g0:g1, :]
    if rows.size == 0:
      continue
    for s0, s1 in sample_slices:
        block = rows[:, s0:s1]
        if block.size == 0:
            continue
        prof = np.nanmean(block, axis=0)   # mean over regions
        v = np.nanmax(prof)
        if np.isfinite(v) and v > maxv:
            maxv = float(v)

print(maxv)
PY
}

pad_ymax () {
  python - "$1" <<'PY'
import sys, math
m=float(sys.argv[1])
y=m*1.05
y=math.ceil(y*1000)/1000.0
print(y)
PY
}

echo "[INFO] INDIR=$INDIR"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] TAG=$TAG"
echo "[INFO] Regions: ${REGIONS[*]}"
echo "[INFO] DMRSETS: ${DMRSETS[*]}"
echo

for REG in "${REGIONS[@]}"; do
  echo "==== REGION: ${REG} ===="

  M1="$(matrix_path "${DMRSETS[0]}" "$REG")"
  M2="$(matrix_path "${DMRSETS[1]}" "$REG")"

  if [[ -z "$M1" || -z "$M2" ]]; then
    echo "[SKIP] missing matrices for region=$REG"
    echo "       expected like: ${INDIR}/C-DMR.${REG}.${TAG}.matrix.gz"
    echo "                      ${INDIR}/CG-DMR.${REG}.${TAG}.matrix.gz"
    echo
    continue
  fi

  RDIR="${OUTDIR}/${REG}"
  mkdir -p "$RDIR"

  max1="$(max_mean_profile_from_matrix "$M1")"
  max2="$(max_mean_profile_from_matrix "$M2")"
  max_all="$(python - <<PY
m1=float("${max1}"); m2=float("${max2}")
print(m1 if m1>=m2 else m2)
PY
)"
  yMax="$(pad_ymax "$max_all")"

  echo "[INFO] maxMean(C-DMR)=${max1}  maxMean(CG-DMR)=${max2}  => shared yMax=${yMax}"

  P1="${RDIR}/${DMRSETS[0]}.${REG}.${TAG}.sharedY.pdf"
  P2="${RDIR}/${DMRSETS[1]}.${REG}.${TAG}.sharedY.pdf"

  # re-plot with shared y axis
  plotProfile -m "$M1" -out "$P1" --perGroup \
    --yMin 0 --yMax "$yMax" \
    --plotTitle "${DMRSETS[0]} ${REG} (${TAG})" \
    --numberOfProcessors "$NCPU"

  plotProfile -m "$M2" -out "$P2" --perGroup \
    --yMin 0 --yMax "$yMax" \
    --plotTitle "${DMRSETS[1]} ${REG} (${TAG})" \
    --numberOfProcessors "$NCPU"

  OUTPDF="${RDIR}/C-vs-CG.${REG}.${TAG}.sharedY.pdf"
  if [[ $have_pdfjam -eq 1 ]]; then
    pdfjam "$P1" "$P2" --nup 2x1 --landscape --outfile "$OUTPDF" >/dev/null 2>&1
    echo "[OK] wrote side-by-side: $OUTPDF"
  elif [[ $have_pdfunite -eq 1 ]]; then
    pdfunite "$P1" "$P2" "$OUTPDF"
    echo "[OK] wrote 2-page PDF: $OUTPDF"
  else
    echo "[WARN] no pdfjam/pdfunite found; leaving separate PDFs:"
    echo "       $P1"
    echo "       $P2"
  fi

  echo
done

echo "[DONE] Outputs under: $OUTDIR"
