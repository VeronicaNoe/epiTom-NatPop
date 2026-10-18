#!/usr/bin/env bash
BASE="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly"
EPI_DIR="${BASE}/ai_DMR-classification/bd_merged-methylation"
GEN_DIR="${BASE}/ak_genetic-variants"
COMP_DIR="${BASE}/al_DMR-gen-var_denovo-pedigree_comparison"
RESULTS_DIR="${BASE}/results"
PLOTS_DIR="${BASE}/plots"
#mkdir -p "${COMP_DIR}" "${RESULTS_DIR}" "${PLOTS_DIR}"

WINDOW_BP="${WINDOW_BP:-1000}"
SNP_VCF="${GEN_DIR}/cervil_SNPs.GQ20.DP10.vcf.gz"
SV_VCF="${GEN_DIR}/cervil_SV.PASS.SUP5.vcf.gz"
DMR_VCFS=(
  "${EPI_DIR}/00.1_dmr_binary_states.C-DMR.vcf.gz"
  "${EPI_DIR}/00.1_dmr_binary_states.CG-DMR.vcf.gz"
)

echo "Output folder: ${COMP_DIR}"
echo "WINDOW_BP=${WINDOW_BP}"

# Genetic variants -> BED
bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\n' "${SNP_VCF}" \
  | awk 'BEGIN{OFS="\t"} {id=$4; if(id==".") id=$1":"$2"-"$3":SNP"; print $1,$2,$3,id,"SNP",$5,$6}' \
  | sort -k1,1V -k2,2n > "${COMP_DIR}/01_snps.bed"

bcftools query -f '%CHROM\t%POS0\t%END\t%ID\t%SVTYPE\n' "${SV_VCF}" \
    | awk 'BEGIN{OFS="\t"} {
        end=$3;
        if(end=="." || end=="" || end < $2+1) end=$2+1;
        id=$4; if(id==".") id=$1":"$2"-"end":SV";
        svtype=$5; if(svtype=="." || svtype=="") svtype="SV";
        print $1,$2,end,id,svtype,".","."
      }' \
    | sort -k1,1V -k2,2n \
    > "${COMP_DIR}/02_svs.bed"

cat "${COMP_DIR}/01_snps.bed" "${COMP_DIR}/02_svs.bed" \
  | sort -k1,1V -k2,2n > "${COMP_DIR}/03_genetic_variants.all.bed"

# DMR synthetic VCF -> 1-bp BED sites
for VCF in "${DMR_VCFS[@]}"; do
  DMR_TYPE=$(basename "${VCF}" | sed 's/00.1_dmr_binary_states.//; s/.vcf.gz//')
  bcftools query -f '%CHROM\t%POS0\t%ID[\t%GT]\n' "${VCF}" \
    | awk -v dmr="${DMR_TYPE}" 'BEGIN{OFS="\t"} {print $1,$2,$2+1,$3,dmr}' \
    | sort -k1,1V -k2,2n > "${COMP_DIR}/04_dmr_sites.${DMR_TYPE}.bed"
done

cat "${COMP_DIR}"/04_dmr_sites.*.bed \
  | sort -k1,1V -k2,2n > "${COMP_DIR}/05_dmr_sites.all.bed"

# Intersections
bedtools intersect \
  -a "${COMP_DIR}/05_dmr_sites.all.bed" \
  -b "${COMP_DIR}/03_genetic_variants.all.bed" \
  -wa -wb -sorted \
  > "${COMP_DIR}/06_dmr_genvar.direct_overlap.tsv"

bedtools window \
  -a "${COMP_DIR}/05_dmr_sites.all.bed" \
  -b "${COMP_DIR}/03_genetic_variants.all.bed" \
  -w "${WINDOW_BP}" \
  > "${COMP_DIR}/07_dmr_genvar.window_${WINDOW_BP}bp.tsv"

