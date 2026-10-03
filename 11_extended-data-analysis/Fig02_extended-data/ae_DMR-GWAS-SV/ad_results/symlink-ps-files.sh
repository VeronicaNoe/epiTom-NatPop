#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

BASE="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results"
LIST="$BASE/all.sig.files.txt"

# ----------------------------
# 1. Delete old sig/*.ps.gz links fast
# ----------------------------
echo "Deleting old sig/*.ps.gz ..."
for dmr in C-DMR CG-DMR; do
  for c in ch{01..12}; do
    sigdir="$BASE/$dmr/$c/sig"
    [ -d "$sigdir" ] || continue
    rm -f "$sigdir"/*.ps.gz
  done
done

# ----------------------------
# 2. Recreate sig/*.ps.gz from current all.sig.files.txt
# ----------------------------
echo "Recreating sig/*.ps.gz from $LIST ..."
sed 's/\r$//' "$LIST" | while IFS= read -r id; do
  [ -n "$id" ] || continue

  chr="${id%%_*}"
  rest="${id#*_}"
  dmr="${rest%%_*}"

  src="$BASE/$dmr/$chr/toBrick/${id}.ps.gz"
  dst="$BASE/$dmr/$chr/sig/${id}.ps.gz"

  if [ -f "$src" ]; then
    ln -sfn "../toBrick/${id}.ps.gz" "$dst"
  else
    echo "Missing: $src"
  fi
done

# ----------------------------
# 3. Verify counts
# ----------------------------
echo
echo "Checking mQTL vs ps.gz counts in sig ..."
for dmr in C-DMR CG-DMR; do
  for c in ch{01..12}; do
    sigdir="$BASE/$dmr/$c/sig"
    [ -d "$sigdir" ] || continue

    n_mqtl=$(find "$sigdir" -maxdepth 1 -type f -name '*.mQTL' | wc -l)
    n_ps=$(find "$sigdir" -maxdepth 1 -name '*.ps.gz' | wc -l)

    status="OK"
    if [ "$n_mqtl" -ne "$n_ps" ]; then
      status="MISMATCH"
    fi

    printf "%-6s %-4s mQTL=%-6s ps.gz=%-6s %s\n" "$dmr" "$c" "$n_mqtl" "$n_ps" "$status"
  done
done
