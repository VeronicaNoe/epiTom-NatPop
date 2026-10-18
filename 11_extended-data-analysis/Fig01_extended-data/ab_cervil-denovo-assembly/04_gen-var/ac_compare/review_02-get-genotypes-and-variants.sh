#!/usr/bin/env bash
PROJ="/work/project/ijpb/vibanez/14_cervil-project/03_fast_reference"
OUT="${PROJ}/ac_compare"

module load bioinfo/Bcftools/1.21
module load bioinfo/tabix/1.18
#
bcftools view \
  -f PASS \
  -Oz \
  -o "${PROJ}/ab_run-sniffles/all-samples_joint.PASS.vcf.gz" \
  "${PROJ}/ab_run-sniffles/all-samples_joint.vcf.gz"

tabix -f -p vcf "${PROJ}/ab_run-sniffles/all-samples_joint.PASS.vcf.gz"

# small variants: merged PASS-only DeepVariant VCF
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  "${OUT}/all-samples.deepvariant.merged.vcf.gz" \
  > "${OUT}/smallvars.gt.tsv"

# SVs: PASS-only joint Sniffles VCF
bcftools query \
  -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%GT]\n' \
  "${PROJ}/ab_run-sniffles/all-samples_joint.PASS.vcf.gz" \
  > "${OUT}/sv.gt.tsv"

python3 review_classify-variants.py small "${OUT}/smallvars.gt.tsv" "${OUT}/smallvars.classified.tsv"
python3 review_classify-variants.py sv    "${OUT}/sv.gt.tsv"       "${OUT}/sv.classified.tsv"
