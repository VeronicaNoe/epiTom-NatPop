#!/usr/bin/env bash
###############################################################################
# Build one closest-TE record for every annotated tomato gene.
#
# IMPORTANT:
# allGenes.bed is currently:
#   chromosome, start, end, strand, geneName
#
# For bedtools closest -D a to use the gene strand correctly, the gene file is
# converted to standard BED6:
#   chromosome, start, end, geneName, score, strand
###############################################################################

GENES="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/allGenes.bed"
TES="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/allTEs.bed"

OUTDIR="/mnt/disk2/vibanez/11_extended-data-analysis/Fig03_extended-data"

GENES_BED6="${OUTDIR}/extended-data-fig03e_all-genes_BED6.bed"
TES_FIXED="${OUTDIR}/extended-data-fig03e_all-TEs_fixed.bed"
GENES_SORTED="${OUTDIR}/extended-data-fig03e_all-genes_BED6_sorted.bed"
TES_SORTED="${OUTDIR}/extended-data-fig03e_all-TEs_sorted.bed"
RAW_CLOSEST="${OUTDIR}/extended-data-fig03e_all-genes_closest-TE_raw.tsv"
FINAL_CLOSEST="${OUTDIR}/extended-data-fig03e_all-genes_closest-TE.tsv"

mkdir -p "${OUTDIR}"

command -v bedtools >/dev/null 2>&1 || {
  echo "ERROR: bedtools is not available in PATH" >&2
  exit 1
}

[[ -s "${GENES}" ]] || {
  echo "ERROR: missing or empty gene BED: ${GENES}" >&2
  exit 1
}

[[ -s "${TES}" ]] || {
  echo "ERROR: missing or empty TE BED: ${TES}" >&2
  exit 1
}


echo "[1/6] Converting genes to standard BED6"

awk -F'\t' 'BEGIN {
  OFS = "\t"
}

NF != 5 {
  print "ERROR: gene line " NR " has " NF " fields: " $0 > "/dev/stderr"
  exit 1
}

$1 == "" ||
$2 !~ /^[0-9]+$/ ||
$3 !~ /^[0-9]+$/ ||
$2 > $3 ||
($4 != "+" && $4 != "-") ||
$5 == "" {
  print "ERROR: invalid gene line " NR ": " $0 > "/dev/stderr"
  exit 1
}

{
  # BED6: chromosome, start, end, name, score, strand
  print $1, $2, $3, $5, ".", $4
}
' "${GENES}" > "${GENES_BED6}"


echo "[2/6] Normalizing TE annotations"

awk -F'\t' 'BEGIN {
  OFS = "\t"
}

NF != 6 {
  print "ERROR: TE line " NR " has " NF " fields: " $0 > "/dev/stderr"
  exit 1
}

$1 == "" ||
$2 !~ /^[0-9]+$/ ||
$3 !~ /^[0-9]+$/ ||
$2 > $3 ||
($4 != "+" && $4 != "-") {
  print "ERROR: invalid TE line " NR ": " $0 > "/dev/stderr"
  exit 1
}

{
  if ($5 == "") $5 = "."
  if ($6 == "") $6 = "."

  print $1, $2, $3, $4, $5, $6
}
' "${TES}" > "${TES_FIXED}"


echo "[3/6] Sorting genes and TEs"

bedtools sort -i "${GENES_BED6}" > "${GENES_SORTED}"
bedtools sort -i "${TES_FIXED}" > "${TES_SORTED}"


echo "[4/6] Finding the closest TE for every gene"

# A has 6 fields, B has 6 fields, and -D adds one field: 13 total.
bedtools closest \
  -a "${GENES_SORTED}" \
  -b "${TES_SORTED}" \
  -D a \
  -t first \
  > "${RAW_CLOSEST}"


echo "[5/6] Adding the output header"

{
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "geneChr" \
    "geneStart" \
    "geneEnd" \
    "geneName" \
    "geneScore" \
    "geneStrand" \
    "teChr" \
    "teStart" \
    "teEnd" \
    "teStrand" \
    "teIdentity" \
    "teSuperfamily" \
    "signedDistance"

  awk -F'\t' 'BEGIN {
    OFS = "\t"
  }

  NF != 13 {
    printf "ERROR: expected 13 columns but found %d at raw line %d\n", NF, NR > "/dev/stderr"
    exit 1
  }

  {
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13
  }
  ' "${RAW_CLOSEST}"
} > "${FINAL_CLOSEST}"


echo "[6/6] Final validation"

N_GENES="$(wc -l < "${GENES_SORTED}")"
N_CLOSEST="$(awk 'END {print NR - 1}' "${FINAL_CLOSEST}")"
N_OVERLAP="$(awk -F'\t' 'NR > 1 && $13 == 0 {n++} END {print n + 0}' "${FINAL_CLOSEST}")"
N_NO_TE="$(awk -F'\t' 'NR > 1 && ($7 == "." || $8 == -1) {n++} END {print n + 0}' "${FINAL_CLOSEST}")"

if [[ "${N_GENES}" -ne "${N_CLOSEST}" ]]; then
  echo "ERROR: gene count (${N_GENES}) differs from closest-TE rows (${N_CLOSEST})" >&2
  exit 1
fi

echo
echo "Completed successfully"
echo "Genes:                         ${N_GENES}"
echo "Closest-TE rows:               ${N_CLOSEST}"
echo "Genes overlapping a TE:        ${N_OVERLAP}"
echo "Genes without a TE chr hit:    ${N_NO_TE}"
echo "Output: ${FINAL_CLOSEST}"
