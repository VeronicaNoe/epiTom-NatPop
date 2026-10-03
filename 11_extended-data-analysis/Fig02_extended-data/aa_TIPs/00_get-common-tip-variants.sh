WDIR="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS"
TIP_RAW="/mnt/disk2/vibanez/otherAnalysis/02_TIPs/TIPs.vcf.gz"
PLINK="/home/vibanez/anaconda3/envs/tls/bin/plink"

OUT_VCF="$WDIR/tables/TIPs.common.ldpruned.vcf.gz"
TMPDIR="$WDIR/tables/_tip_prune_tmp"
mkdir -p "$WDIR/tables" "$TMPDIR"

# 1) Build a PLINK-readable TIP VCF (sanitize REF/ALT only; keep IDs and GTs)
TIP_PLINK_VCF="$TMPDIR/TIPs.plink_sane.vcf.gz"
zcat "$TIP_RAW" \
  | awk 'BEGIN{FS=OFS="\t"}
         /^##/ {print; next}
         /^#CHROM/ {print; next}
         { $4="A"; $5="T"; print }' \
  | bgzip -c > "$TIP_PLINK_VCF"
tabix -f -p vcf "$TIP_PLINK_VCF"

# 2) Import to PLINK
"$PLINK" \
  --vcf "$TIP_PLINK_VCF" \
  --double-id \
  --allow-extra-chr \
  --biallelic-only strict \
  --make-bed \
  --out "$TMPDIR/tip"

# 3) Keep common variants + low missingness
"$PLINK" \
  --bfile "$TMPDIR/tip" \
  --allow-extra-chr \
  --maf 0.05 \
  --geno 0.1 \
  --make-bed \
  --out "$TMPDIR/tip.common"

# 4) LD pruning
"$PLINK" \
  --bfile "$TMPDIR/tip.common" \
  --allow-extra-chr \
  --indep-pairwise 1000 50 0.2 \
  --out "$TMPDIR/tip.common"

# 5) Filter ORIGINAL TIP VCF by pruned IDs (keeps original REF/ALT strings)
KEEP_IDS="$TMPDIR/tip.common.prune.in"
zcat "$TIP_RAW" \
  | awk 'BEGIN{FS=OFS="\t"}
         NR==FNR {keep[$1]=1; next}
         /^#/ {print; next}
         ($3 in keep) {print}' \
      "$KEEP_IDS" - \
  | bgzip -c > "$OUT_VCF"
tabix -f -p vcf "$OUT_VCF"

echo "Wrote: $OUT_VCF"
echo "Counts:"
echo -n "  raw TIP records: "; zcat "$TIP_RAW" | grep -vc '^#'
echo -n "  pruned TIP records: "; zcat "$OUT_VCF" | grep -vc '^#'
