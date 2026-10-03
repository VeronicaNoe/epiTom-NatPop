for file in *CG-DMR_sigSNPs.mQTL; do
	awk '{OFS="\t"} NR>1 {print $1,$4,$5,$7,$6}' "$file" | sed 's/:/\t/g' >> CG-DMR_sigSNPs.merged
done

for file in *C-DMR_sigSNPs.mQTL; do
         awk '{OFS="\t"} NR>1 {print $1,$4,$5,$7,$6}' "$file" | sed 's/:/\t/g' >> C-DMR_sigSNPs.merged
done
##
# 1) Extract SV ID -> CHROM POS from the VCF
#bcftools query -f '%ID\t%CHROM\t%POS\n' /mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/aa_markers/graphpanSV.sharedSamples.vcf.gz \
#  | sort -k1,1 > sv.id2coord.tsv

# 2) Sort your merged files by SV ID
sort -k1,1 CG-DMR_sigSNPs.merged > CG-DMR_sigSNPs.sorted
sort -k1,1 C-DMR_sigSNPs.merged  > C-DMR_sigSNPs.sorted

# 3) Join by SV ID
join -t $'\t' -1 1 -2 1 sv.id2coord.tsv CG-DMR_sigSNPs.sorted > CG-DMR_sigSNPs.withcoord
join -t $'\t' -1 1 -2 1 sv.id2coord.tsv C-DMR_sigSNPs.sorted  > C-DMR_sigSNPs.withcoord


awk 'BEGIN{OFS="\t"} {printf "%02d\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $3+1, $1, $4, "ch"$5"_"$6"_"$7}' CG-DMR_sigSNPs.withcoord \
  | sortBed -i - \
  | mergeBed -i - -c 4,5,6,6 -o distinct,min,count_distinct,distinct \
  > ../bc_annotation/CG-DMR_sigSNPs.bed

awk 'BEGIN{OFS="\t"} {printf "%02d\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $3+1, $1, $4, "ch"$5"_"$6"_"$7}' C-DMR_sigSNPs.withcoord \
  | sortBed -i - \
  | mergeBed -i - -c 4,5,6,6 -o distinct,min,count_distinct,distinct \
  > ../bc_annotation/C-DMR_sigSNPs.bed
