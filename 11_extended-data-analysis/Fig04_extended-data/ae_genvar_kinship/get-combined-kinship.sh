#!/usr/bin/env bash
# =========================
# INPUTS
# =========================
TFAM="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/DMR_general_leaf_LD.tfam"

SNPs="/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/ab_vcf-metabolome/SNP_general_leaf-metabolome_LD.vcf.gz"
SV="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/aa_markers/graphpanSV.sharedSamples.vcf.gz"
SV_COORD="/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/aa_DMR-GWAS-SV_SL25/sv_id2coord_sl25_primary.tsv"

TIP_RAW="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS/tables/TIPs.common.ldpruned.vcf.gz"
TIP="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS/tables/TIPs.common.ldpruned.fixed.vcf.gz"

OUTDIR="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS/ae_genvar_kinship"
mkdir -p "$OUTDIR"

PLINK="/home/vibanez/anaconda3/envs/tls/bin/plink"
EMMAX_KIN="$HOME/bin/tools/EMMAX/emmax-kin-intel64"

SV_ID_COL=1
SV_CHR_COL=4
SV_POS_COL=5

# =========================
# 0) FIX TIP VCF HEADER
# =========================
zcat "$TIP_RAW" | awk '
BEGIN {
  print "##fileformat=VCFv4.2"
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
}
{ print }
' | bgzip -c > "$TIP"

tabix -f -p vcf "$TIP"

# =========================
# 1) SAMPLE LIST FROM DMR TFAM
# =========================
awk '{print $2}' "$TFAM" > "$OUTDIR/tfam.samples"
awk '{print $1, $2}' "$TFAM" > "$OUTDIR/tfam.keep"

sort "$OUTDIR/tfam.samples" > "$OUTDIR/tfam.samples.sorted"

echo "TFAM samples:"
wc -l "$OUTDIR/tfam.samples"

echo "VCF sample counts:"
bcftools query -l "$SNPs" | wc -l
bcftools query -l "$TIP"  | wc -l
bcftools query -l "$SV"   | wc -l

bcftools query -l "$SNPs" | sort > "$OUTDIR/SNP.samples"
bcftools query -l "$TIP"  | sort > "$OUTDIR/TIP.samples"
bcftools query -l "$SV"   | sort > "$OUTDIR/SV.samples"

echo "Samples in TFAM missing from SNP:"
comm -23 "$OUTDIR/tfam.samples.sorted" "$OUTDIR/SNP.samples" || true

echo "Samples in TFAM missing from TIP:"
comm -23 "$OUTDIR/tfam.samples.sorted" "$OUTDIR/TIP.samples" || true

echo "Samples in TFAM missing from SV:"
comm -23 "$OUTDIR/tfam.samples.sorted" "$OUTDIR/SV.samples" || true

# =========================
# 2) CONVERT EACH VCF TO PLINK
# =========================
"$PLINK" \
  --vcf "$SNPs" \
  --double-id \
  --allow-extra-chr \
  --vcf-half-call missing \
  --make-bed \
  --out "$OUTDIR/SNP_raw"

"$PLINK" \
  --vcf "$TIP" \
  --double-id \
  --allow-extra-chr \
  --vcf-half-call missing \
  --make-bed \
  --out "$OUTDIR/TIP_raw"

"$PLINK" \
  --vcf "$SV" \
  --double-id \
  --allow-extra-chr \
  --vcf-half-call missing \
  --make-bed \
  --out "$OUTDIR/SV_raw"

# =========================
# 3) PREFIX SNP/TIP MARKER IDS
# =========================
cp "$OUTDIR/SNP_raw.bed" "$OUTDIR/SNP.bed"
cp "$OUTDIR/SNP_raw.fam" "$OUTDIR/SNP.fam"

awk 'BEGIN{OFS="\t"} {
  $2 = "SNP_" $1 "_" $4 "_" NR
  print
}' "$OUTDIR/SNP_raw.bim" > "$OUTDIR/SNP.bim"

cp "$OUTDIR/TIP_raw.bed" "$OUTDIR/TIP.bed"
cp "$OUTDIR/TIP_raw.fam" "$OUTDIR/TIP.fam"

awk 'BEGIN{OFS="\t"} {
  $2 = "TIP_" $1 "_" $4 "_" NR
  print
}' "$OUTDIR/TIP_raw.bim" > "$OUTDIR/TIP.bim"

# =========================
# 4) REPLACE SV COORDS WITH SL2.5 COORDS WHEN AVAILABLE
# =========================
cp "$OUTDIR/SV_raw.bed" "$OUTDIR/SV.bed"
cp "$OUTDIR/SV_raw.fam" "$OUTDIR/SV.fam"

rm -f "$OUTDIR/SV_ids_missing_SL25_coord.txt"

awk -v idc="$SV_ID_COL" -v chc="$SV_CHR_COL" -v pc="$SV_POS_COL" '
BEGIN{FS=OFS="\t"}
NR==FNR {
  if ($pc !~ /^[0-9]+$/) next
  chr[$idc]=$chc
  pos[$idc]=$pc
  next
}
{
  oldid=$2
  if (oldid in chr) {
    $1=chr[oldid]
    $4=pos[oldid]
  } else {
    print oldid >> "'"$OUTDIR"'/SV_ids_missing_SL25_coord.txt"
  }
  $2 = "SV_" oldid
  print
}
' "$SV_COORD" "$OUTDIR/SV_raw.bim" > "$OUTDIR/SV.bim"

# =========================
# 5) MERGE SNP + TIP + SV
# =========================
cat > "$OUTDIR/merge_list.txt" <<EOF
$OUTDIR/TIP.bed $OUTDIR/TIP.bim $OUTDIR/TIP.fam
$OUTDIR/SV.bed  $OUTDIR/SV.bim  $OUTDIR/SV.fam
EOF

"$PLINK" \
  --bfile "$OUTDIR/SNP" \
  --merge-list "$OUTDIR/merge_list.txt" \
  --allow-extra-chr \
  --make-bed \
  --out "$OUTDIR/genvar_raw"

# =========================
# 6) KEEP DMR SAMPLES AND FORCE DMR SAMPLE ORDER
# =========================
"$PLINK" \
  --bfile "$OUTDIR/genvar_raw" \
  --keep "$OUTDIR/tfam.keep" \
  --indiv-sort file "$OUTDIR/tfam.keep" \
  --allow-extra-chr \
  --make-bed \
  --out "$OUTDIR/genvar_leaf"

# =========================
# 7) CHECK SAMPLE ORDER
# =========================
awk '{print $1,$2}' "$TFAM" > "$OUTDIR/dmr.order"
awk '{print $1,$2}' "$OUTDIR/genvar_leaf.fam" > "$OUTDIR/genvar.order"

echo "Checking sample order after merge:"
if diff -q "$OUTDIR/dmr.order" "$OUTDIR/genvar.order" >/dev/null; then
  echo "OK: genvar_leaf.fam matches DMR TFAM order"
else
  echo "ERROR: sample order still differs after --indiv-sort"
  echo "Check with:"
  echo "diff -u $OUTDIR/dmr.order $OUTDIR/genvar.order | head -40"
  exit 1
fi

echo "Marker counts:"
wc -l "$OUTDIR/SNP.bim"
wc -l "$OUTDIR/TIP.bim"
wc -l "$OUTDIR/SV.bim"
wc -l "$OUTDIR/genvar_leaf.bim"

# =========================
# 8) MISSINGNESS CHECK
# =========================
"$PLINK" \
  --bfile "$OUTDIR/genvar_leaf" \
  --allow-extra-chr \
  --missing \
  --out "$OUTDIR/genvar_leaf_missingness"

# =========================
# 9) LD PRUNING
# =========================
"$PLINK" \
  --bfile "$OUTDIR/genvar_leaf" \
  --allow-extra-chr \
  --geno 0.2 \
  --maf 0.05 \
  --indep-pairwise 50 5 0.2 \
  --out "$OUTDIR/genvar_leaf_LD"

"$PLINK" \
  --bfile "$OUTDIR/genvar_leaf" \
  --extract "$OUTDIR/genvar_leaf_LD.prune.in" \
  --indiv-sort file "$OUTDIR/tfam.keep" \
  --allow-extra-chr \
  --make-bed \
  --out "$OUTDIR/genvar_leaf_LD"

# =========================
# 10) CHECK FINAL LD SAMPLE ORDER
# =========================
awk '{print $1,$2}' "$OUTDIR/genvar_leaf_LD.fam" > "$OUTDIR/genvar_LD.order"

echo "Checking sample order after LD pruning:"
if diff -q "$OUTDIR/dmr.order" "$OUTDIR/genvar_LD.order" >/dev/null; then
  echo "OK: genvar_leaf_LD.fam matches DMR TFAM order"
else
  echo "ERROR: LD-pruned sample order differs from DMR TFAM"
  echo "Check with:"
  echo "diff -u $OUTDIR/dmr.order $OUTDIR/genvar_LD.order | head -40"
  exit 1
fi

# =========================
# 11) RECODE TO NUMERIC TPED/TFAM
# =========================
"$PLINK" \
  --bfile "$OUTDIR/genvar_leaf_LD" \
  --allow-extra-chr \
  --recode 12 transpose \
  --out "$OUTDIR/genvar_leaf_LD"

# Check TPED/TFAM order one last time
awk '{print $1,$2}' "$OUTDIR/genvar_leaf_LD.tfam" > "$OUTDIR/genvar_LD_tped.order"

echo "Checking final TPED/TFAM sample order:"
if diff -q "$OUTDIR/dmr.order" "$OUTDIR/genvar_LD_tped.order" >/dev/null; then
  echo "OK: genvar_leaf_LD.tfam matches DMR TFAM order"
else
  echo "ERROR: final TPED/TFAM order differs from DMR TFAM"
  exit 1
fi

# =========================
# 12) BUILD EMMAX KINSHIP
# =========================
cd "$OUTDIR"

"$EMMAX_KIN" -v -d 10 genvar_leaf_LD 2>&1 | tee genvar_leaf_LD.BN.emmax-kin.log

"$EMMAX_KIN" -v -s -d 10 genvar_leaf_LD 2>&1 | tee genvar_leaf_LD.IBS.emmax-kin.log

echo "Done."
echo "Final files:"
ls -lh "$OUTDIR"/genvar_leaf_LD*.kinf
