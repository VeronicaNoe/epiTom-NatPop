#!/usr/bin/env bash
# ============================================================
# Paths
# ============================================================
OUTDIR="/mnt/disk2/vibanez/otherAnalysis/12_getting-graphpanSV-sample-names"
GRAPH_DIR="/mnt/disk2/vibanez/otherAnalysis/00_public-data/genotypes/genotypes"
MY_VCF="/mnt/disk2/vibanez/otherAnalysis/08_liftover-SNPs-SL25-to-SL5/SNP_LD.SL5.fix.vcf.gz"
REF="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SL5.0.fasta"
# ============================================================
# Parallelization
# IMPORTANT: N_JOBS * THREADS_PER_JOB should fit your machine
# ============================================================
N_JOBS=2
THREADS_PER_JOB=4
THREADS_MAIN=12
# ============================================================
# Optional chromosome rename maps
# Leave empty ("") if not needed
#
# Format: old_name <tab> new_name
# Example if graph VCF uses chr1 but REF/my VCF use SL5.0ch01:
# chr1    SL5.0ch01
# chr2    SL5.0ch02
# ...
# ============================================================
GRAPH_RENAME_MAP=""
MY_RENAME_MAP=""

# ============================================================
# Output structure
# ============================================================
mkdir -p \
  "$OUTDIR"/00_graph-snps \
  "$OUTDIR"/01_norm \
  "$OUTDIR"/02_shared \
  "$OUTDIR"/03_merged \
  "$OUTDIR"/04_plink \
  "$OUTDIR"/logs

# ============================================================
# Sanity checks
# ============================================================
# conda activate tls
for x in bcftools plink bgzip tabix; do
  command -v "$x" >/dev/null 2>&1 || { echo "[ERROR] Missing: $x";  }
done
[[ -f "$REF" ]] || { echo "[ERROR] REF not found: $REF"; }
[[ -f "$MY_VCF" ]] || { echo "[ERROR] MY_VCF not found: $MY_VCF";  }
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"
[[ -f "${MY_VCF}.tbi" || -f "${MY_VCF}.csi" ]] || bcftools index -t "$MY_VCF"

# ============================================================
# Chromosome list in correct order
# ============================================================
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12)

# ============================================================
# Step 1: extract SNPs from graphpan chr1-12 in parallel
# Keeps only biallelic SNPs
# ============================================================
extract_graph-snps_one() {
  local chr="$1"
  local in_vcf="${GRAPH_DIR}/${chr}.vcf.gz"
  local out_vcf="${OUTDIR}/00_graph-snps/${chr}.snps.vcf.gz"
  [[ -f "$in_vcf" ]] || { echo "[ERROR] Missing $in_vcf"; }
  bcftools view \
    -v snps -m2 -M2 \
    --threads "$THREADS_PER_JOB" \
    -Oz -o "$out_vcf" \
    "$in_vcf"
  bcftools index -t "$out_vcf"
}
export -f extract_graph-snps_one
export GRAPH_DIR OUTDIR THREADS_PER_JOB

printf "%s\n" "${CHRS[@]}" | xargs -I{} -n1 -P "$N_JOBS" bash -c 'extract_graph-snps_one "$@"' _ {}

# ============================================================
# Step 2: concat graphpan SNP VCFs
# NOTE: concat, not merge, because these are chromosome-split files
# ============================================================
GRAPH_SNP_RAW="${OUTDIR}/00_graph-snps/graphpan.snps.raw.vcf.gz"

bcftools concat \
  --threads "$THREADS_MAIN" \
  -Oz -o "$GRAPH_SNP_RAW" \
  "${OUTDIR}/00_graph-snps/chr1.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr2.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr3.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr4.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr5.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr6.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr7.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr8.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr9.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr10.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr11.snps.vcf.gz" \
  "${OUTDIR}/00_graph-snps/chr12.snps.vcf.gz"

bcftools index -t "$GRAPH_SNP_RAW"

# ============================================================
# Helper: normalize VCF
# - optional chr rename
# - left-normalize vs REF
# - remove REF mismatches (-c x)
# - split multiallelics if any
# - keep only biallelic SNPs
# - force ID = CHROM:POS:REF:ALT
# ============================================================
normalize_vcf() {
  local in_vcf="$1"
  local out_vcf="$2"
  local rename_map="$3"
  if [[ -n "$rename_map" ]]; then
    bcftools annotate --rename-chrs "$rename_map" -Ou "$in_vcf" \
    | bcftools norm -f "$REF" -c x -m -any -Ou \
    | bcftools view -v snps -m2 -M2 -Ou \
    | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads "$THREADS_MAIN" -Oz -o "$out_vcf"
  else
    bcftools norm -f "$REF" -c x -m -any -Ou "$in_vcf" \
    | bcftools view -v snps -m2 -M2 -Ou \
    | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads "$THREADS_MAIN" -Oz -o "$out_vcf"
  fi
  bcftools index -t "$out_vcf"
}

export REF THREADS_MAIN
export -f normalize_vcf

# ============================================================
# Step 3: normalize graphpan SNP VCF
# ============================================================
GRAPH_SNP_NORM="${OUTDIR}/01_norm/graphpan.snps.norm.vcf.gz"
normalize_vcf "$GRAPH_SNP_RAW" "$GRAPH_SNP_NORM" "$GRAPH_RENAME_MAP"

# ============================================================
# Step 4: normalize your named SNP VCF
# ============================================================
MY_SNP_NORM="${OUTDIR}/01_norm/my.named.snps.norm.vcf.gz"
normalize_vcf "$MY_VCF" "$MY_SNP_NORM" "$MY_RENAME_MAP"
##
OUTDIR="/mnt/disk2/vibanez/otherAnalysis/12_getting-graphpanSV-sample-names"
MY_VCF="/mnt/disk2/vibanez/otherAnalysis/08_liftover-SNPs-SL25-to-SL5/SNP_LD.SL5.fix.vcf.gz"
REF="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SL5.0.fasta"
mkdir -p "$OUTDIR/01_norm"
MY_GT_ONLY="${OUTDIR}/01_norm/my.named.GTonly.vcf.gz"
MY_SNP_NORM="${OUTDIR}/01_norm/my.named.snps.norm.vcf.gz"
zcat "$MY_VCF" \
| awk 'BEGIN{FS=OFS="\t"}
  /^#/ {print; next}
  {
    $9="GT";
    for(i=10;i<=NF;i++){
      split($i,a,":");
      $i=a[1];
    }
    print
  }' \
| bgzip -@ 12 -c > "$MY_GT_ONLY"

tabix -p vcf "$MY_GT_ONLY"
bcftools norm -f "$REF" -c x -m -any -Ou "$MY_GT_ONLY" \
| bcftools view -v snps -m2 -M2 -Ou \
| bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 12 -Oz -o "$MY_SNP_NORM"

bcftools index -t "$MY_SNP_NORM"
##
# ============================================================
# Step 5: get exact shared SNP IDs
# Since IDs are CHROM:POS:REF:ALT, this is exact matching
# ============================================================
GRAPH_IDS="${OUTDIR}/02_shared/graphpan.ids.txt"
MY_IDS="${OUTDIR}/02_shared/my.ids.txt"
SHARED_IDS="${OUTDIR}/02_shared/shared.ids.txt"

bcftools query -f '%ID\n' "$GRAPH_SNP_NORM" | LC_ALL=C sort -u > "$GRAPH_IDS"
bcftools query -f '%ID\n' "$MY_SNP_NORM"    | LC_ALL=C sort -u > "$MY_IDS"
LC_ALL=C comm -12 "$GRAPH_IDS" "$MY_IDS" > "$SHARED_IDS"

echo "[INFO] graphpan normalized SNPs: $(wc -l < "$GRAPH_IDS")"
echo "[INFO] my normalized SNPs:       $(wc -l < "$MY_IDS")"
echo "[INFO] exact shared SNPs:        $(wc -l < "$SHARED_IDS")"

# ============================================================
# Step 6: restrict both VCFs to the exact same shared SNPs
# ============================================================
GRAPH_SHARED="${OUTDIR}/02_shared/graphpan.shared.vcf.gz"
MY_SHARED="${OUTDIR}/02_shared/my.named.shared.vcf.gz"

bcftools view \
  --threads "$THREADS_MAIN" \
  -i "ID=@${SHARED_IDS}" \
  -Oz -o "$GRAPH_SHARED" \
  "$GRAPH_SNP_NORM"
bcftools index -t "$GRAPH_SHARED"

bcftools view \
  --threads "$THREADS_MAIN" \
  -i "ID=@${SHARED_IDS}" \
  -Oz -o "$MY_SHARED" \
  "$MY_SNP_NORM"
bcftools index -t "$MY_SHARED"

echo "[INFO] graphpan shared SNP records: $(bcftools view -H "$GRAPH_SHARED" | wc -l)"
echo "[INFO] my shared SNP records:       $(bcftools view -H "$MY_SHARED" | wc -l)"

# ============================================================
# Step 7: merge both shared VCFs
# IMPORTANT:
# - if sample names overlap between files, bcftools merge will fail
# - if that happens, reheader one dataset with prefixes first
# ============================================================
MERGED_VCF="${OUTDIR}/03_merged/graphpan_plus_my.shared.vcf.gz"

bcftools merge \
  --threads "$THREADS_MAIN" \
  -m none \
  -Oz -o "$MERGED_VCF" \
  "$GRAPH_SHARED" \
  "$MY_SHARED"

bcftools index -t "$MERGED_VCF"

# ============================================================
# Step 8: convert to PLINK BED
# ============================================================
PLINK_PREFIX="${OUTDIR}/04_plink/graphpan_plus_my.shared"

plink \
  --vcf "$MERGED_VCF" \
  --vcf-half-call missing \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --threads "$THREADS_MAIN" \
  --out "$PLINK_PREFIX"
  
# ============================================================
# Step 9: compute 1-IBS distance matrix
# Output:
# - .mdist
# - .mdist.id
# ============================================================
plink \
  --bfile "$PLINK_PREFIX" \
  --allow-extra-chr \
  --distance square 1-ibs flat-missing \
  --threads "$THREADS_MAIN" \
  --out "$PLINK_PREFIX"

# ============================================================
# Optional: sample lists for later matching
# ============================================================
bcftools query -l "$GRAPH_SHARED" > "${OUTDIR}/05_matches/graphpan.samples.txt"
bcftools query -l "$MY_SHARED"    > "${OUTDIR}/05_matches/my.samples.txt"
bcftools query -l "$MERGED_VCF"   > "${OUTDIR}/05_matches/merged.samples.txt"

echo "[DONE]"
echo "[OUT] Shared IDs:        $SHARED_IDS"
echo "[OUT] Graph shared VCF:  $GRAPH_SHARED"
echo "[OUT] My shared VCF:     $MY_SHARED"
echo "[OUT] Merged VCF:        $MERGED_VCF"
echo "[OUT] PLINK BED:         ${PLINK_PREFIX}.bed/.bim/.fam"
echo "[OUT] 1-IBS distance:    ${PLINK_PREFIX}.mdist + ${PLINK_PREFIX}.mdist.id"
