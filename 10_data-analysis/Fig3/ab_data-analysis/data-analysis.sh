####
## merge H2 from all the results
Rscript --vanilla ~/bin/fig3_00.1_merge-h2.R
## do the pltos
Rscript --vanilla ~/bin/fig3_00.2_merge-GIF-MAf-H2.R
## from the snp block, get cis and trans annotation
Rscript --vanilla ~/bin/fig3_01.0_get-cis-trans-association.R
## combine pvalues
python ~/bin/fig3_02.1_combine_GWAS-pvalues-chisq.py
# get the manhattan with combined pvalues
python ~/bin/fig3_02.2_get_manhattan-combinedChisq.py
# get the heatmap of pvalues
# HERE we need to run the script in the dir with the files to use
# apparently I splic CG and C sig to do each
python ~/bin/fig3_02.3_get_SNPs-DMRs-pvalues-heatmap.py
### get most sig SNPs
python ~/bin/fig3_02.4_get_top_1_SNPs.py
### get the genes closest to those SNPs
ANNODIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
INDIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results"
cat ${INDIR}/02.2_CG-DMR_top_1_percent.csv | sed 's/:/\t/g' |\
sed 's/,/\t/g'|awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $2+1, $3, $4, $5}'  |\
intersectBed -a - -b ${ANNODIR}/allGeneNames-Function.bed -wo > ${INDIR}/02.3_CG-DMR.closest-gene_top-1-SNPs.bed

cat ${INDIR}/02.2_C-DMR_top_1_percent.csv | sed 's/:/\t/g' |\
sed 's/,/\t/g'|awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $2+1, $3, $4, $5}'  |\
intersectBed -a - -b ${ANNODIR}/allGeneNames-Function.bed -wo > ${INDIR}/02.3_C-DMR.closest-gene_top-1-SNPs
##
Rscript --vanilla ~/bin/fig3_03.1_get_SNPs_over_epiGenes.R
## split each snp over epigene
# outpit in associated-SNP-over-epigenes
# old dir dmrAssSNP...
# we have other dir in  /mnt/disk3/vibanez/DMR-GWAS/results/snpOverEpiGene
# the gene list seems larger, ah, maybe is because we started with top 5% of ass snps
bash ~/bin/check_top5chisq.sh

### get the meth levels of the DMRs
# in this script we merge data from SNPs, the meth levels and status and the snpEffect over the gene
Rscript --vanilla ~/bin/fig3_04.1_get_DMRs_linked2SNPs_overEpiGenes.R
# also, we annotate the SNP
# 03.3_annottate_SNPs.R
# snpEff dir:
Mcclintock:/mnt/data6/vibanez/SNPs/vcfiles/snpEff
macos:/Users/macbook/Desktop/tmo/Fig3/data/snpEff
# I need to clean this. we need to process with snpEffect before R plot and
# I'm not sure why we are annotating the snps again.
# this file come from the python scrip 04.0
#'03.9_SNPs_over_epiGenes_annotated.tsv')

# 04.4_candidate_genes_dmr_meth (2 files)
I don't understand why we use KO.
Also, use the ko data from other figure, don't duplicate files
