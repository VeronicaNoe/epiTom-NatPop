###
# allele frequency spectrum
Rscript --vanilla ~/bin/fig2_01.1_allele-frequency-spectrum.R
###
# LD decay
Rscript --vanulla ~/bin/fig2_02.1_LD-decay.R
###
# (epi)genetic diversity
# the plots are with 20 accessiones from each group
# this was because the value could change if they have unbalance number of individuals
Rscript --vanilla ~/bin/fig2_03.1_epi-PI.R
####
# pairwise differences
#https://www.cog-genomics.org/plink/2.0/diff#sample_diff
  ~/bin/pairwise-divergence.sh
  ~/bin/pairwise-divergence_DMRs.sh
    plink2 --vcf  --make-bed --allow-extra-chr --out test
    plink2 --bfile test --double-id --allow-extra-chr --threads 10 --sample-diff pairwise ids=P12_leaf  P13_leaf  P1_leaf --out test

   # output
   #FID1   IID1    FID2    IID2    OBS_CT  DIFF_CT
   #FID1	maybefid, fid	FID of first sample in current pair
   #IID1	(required)	IID of first sample in current pair
   #FID2	maybefid, fid	FID of second sample in current pair
   #IID2	(required)	IID of second sample in current pair
   #OBS_CT	nobs	Number of genotype/dosage pairs considered
   #DIFF_CT	diff	# of genotype/dosage discordances

