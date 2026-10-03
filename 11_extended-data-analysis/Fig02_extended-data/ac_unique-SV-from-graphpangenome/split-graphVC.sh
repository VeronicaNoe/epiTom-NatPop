WD="/mnt/disk2/vibanez/otherAnalysis/08_liftover-SNPs-SL25-to-SL5"
V_GONLY="$WD/graph_only_all_not_in_my.fix.SL5.vcf.gz"

OUTROOT="/mnt/disk2/vibanez/otherAnalysis/10_unique-SV-from-graphpangenome"
UNIQ_VCF_DIR="$OUTROOT/aa_graphpan-unique-vcf"
mkdir -p "$UNIQ_VCF_DIR"

for c in 0 1 2 3 4 5 6 7 8 9 10 11 12; do
  bcftools view --threads 24 -r "$c" -Oz -o "$UNIQ_VCF_DIR/chr${c}.vcf.gz" "$V_GONLY"
  bcftools index --threads 24 -t "$UNIQ_VCF_DIR/chr${c}.vcf.gz"
done
