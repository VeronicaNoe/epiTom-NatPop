#!/usr/bin/env bash
INDIR="/mnt/disk2/vibanez/otherAnalysis/12_getting-graphpanSV-sample-names/11_sv-final"
OUTDIR="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/aa_markers"
GP_SV="${INDIR}/graphpanSV.selected.renamed.vcf.gz"
SNPs="/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/ba_markers/SNPs.LD.renamed.vcf.gz"
SNP_for_SV="$OUTDIR/SNPs.sharedSamples-with-SV.vcf.gz"
SV_shared="$OUTDIR/graphpanSV.sharedSamples.vcf.gz"
SV_PREFIX="$OUTDIR/graphpanSV"
SV_LD_PREFIX="$OUTDIR/graphpanSV.LD"
SNP_PREFIX="$OUTDIR/SNPs.sharedSamples-with-SV"
SNP_LD_PREFIX="$OUTDIR/SNPs.sharedSamples-with-SV.LD"

THREADS=40

# ------------------------------------------------------------
# 1) build shared sample list
# ------------------------------------------------------------
echo "[ STEP 1 ]: build shared sample list"
bcftools query -l "$GP_SV" > "$OUTDIR/sv.samples.inorder.txt"
bcftools query -l "$SNPs"  | sort -u > "$OUTDIR/snp.samples.sorted.txt"

awk 'NR==FNR{a[$1]=1; next} ($1 in a)' \
  "$OUTDIR/snp.samples.sorted.txt" \
  "$OUTDIR/sv.samples.inorder.txt" \
  > "$OUTDIR/shared.samples.raw"

# remove samples you do not want
grep -vxE 'TS-244|S_pimLA1578_1' "$OUTDIR/shared.samples.raw" > "$OUTDIR/shared.samples"

echo "[INFO] Number of shared samples:"
wc -l "$OUTDIR/shared.samples"

#echo "[INFO] Shared samples:"
#cat "$OUTDIR/shared.samples"

# ------------------------------------------------------------
# 2) subset both VCFs to the same samples
# ------------------------------------------------------------
echo "[ STEP 2 ]: subset both VCFs to the same samples"
bcftools view \
  -S "$OUTDIR/shared.samples" \
  -Oz -o "$SV_shared" \
  "$GP_SV"

bcftools index -t "$SV_shared"

bcftools view \
  -S "$OUTDIR/shared.samples" \
  -Oz -o "$SNP_for_SV" \
  "$SNPs"

bcftools index -t "$SNP_for_SV"

# ------------------------------------------------------------
# 3) prepare SV markers for EMMAX GWAS
# ------------------------------------------------------------
echo "[ STEP 3 ]: prepare SV markers for EMMAX GWAS"
plink \
  --vcf "$SV_shared" \
  --vcf-half-call missing \
  --threads "$THREADS" \
  --double-id \
  --allow-extra-chr \
  --maf 0.05 \
  --make-bed \
  --out "$SV_PREFIX"

plink \
  --bfile "$SV_PREFIX" \
  --allow-extra-chr \
  --indep-pairwise 100 100 0.2 \
  --out "$SV_PREFIX"

plink \
  --bfile "$SV_PREFIX" \
  --extract "${SV_PREFIX}.prune.in" \
  --make-bed \
  --allow-extra-chr \
  --out "$SV_LD_PREFIX"

# EMMAX marker input
plink \
  --bfile "$SV_LD_PREFIX" \
  --recode12 \
  --threads "$THREADS" \
  --double-id \
  --allow-extra-chr \
  --output-missing-genotype 0 \
  --transpose \
  --out "$SV_LD_PREFIX"

# ------------------------------------------------------------
# 4) prepare SNPs for kinship
# ------------------------------------------------------------
echo "[ STEP 4 ]: prepare SNPs for kinship"
plink \
  --vcf "$SNP_for_SV" \
  --vcf-half-call missing \
  --threads "$THREADS" \
  --double-id \
  --allow-extra-chr \
  --maf 0.05 \
  --make-bed \
  --out "$SNP_PREFIX"

# SNP LD pruning for kinship
plink \
  --bfile "$SNP_PREFIX" \
  --allow-extra-chr \
  --indep-pairwise 100 100 0.2 \
  --out "$SNP_PREFIX"

plink \
  --bfile "$SNP_PREFIX" \
  --extract "${SNP_PREFIX}.prune.in" \
  --make-bed \
  --allow-extra-chr \
  --out "$SNP_LD_PREFIX"

plink --bfile "$SNP_LD_PREFIX" --recode12 --threads 40 --double-id --allow-extra-chr \
  --output-missing-genotype 0 --transpose --out $SNP_LD_PREFIX

# ------------------------------------------------------------
# 5) kinship from SNPs
# ------------------------------------------------------------
echo "[ STEP 5 ]:  kinship from SNPs"
~/bin/EMMAX/emmax-kin-intel64 -v -d 10 "$SNP_LD_PREFIX"     # BN kinship
~/bin/EMMAX/emmax-kin-intel64 -v -s -d 10 "$SNP_LD_PREFIX"  # IBS kinship

echo "[DONE]"
echo "[SV markers for EMMAX]   ${SV_LD_PREFIX}.tped / ${SV_LD_PREFIX}.tfam"
echo "[SNP kinship prefix]     ${SNP_LD_PREFIX}"
