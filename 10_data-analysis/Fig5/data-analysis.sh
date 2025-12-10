#!/bin/bash
## Get number of metabolites with QTLs per molecular marker and kinship
## PLot the heritability of each marker

Rscript --vanilla ~/bin/fig5_01.0_general-description.R

## Merge methylation status, methylation levels, metabolite levels
Rscript --vanilla ~/bin/fig5_01.2_merge_methLevels.R
Rscript --vanilla ~/bin/fig5_01.3_merge_meth-levels-status-metaLevels.R

## Get correlations
Rscript --vanilla ~/bin/fig5_01.4_methylation-metabolite-correlation.R
## get the closest gene of the QTL and the QTLs over TEs and the closest gene
bash ~/bin/fig5_01.5_get-closest-gene.sh
Rscript --vanilla ~/bin/fig5_01.6_merge-closest-gene.R

# get gene expression for natural populaions and the ddm1 ko
# merge RNA counts ( split 05.0 )
Rscript --vanilla ~/bin/fig5_01.7_merge_leaf-transcriptome.R
# merge all info ( from te 05.0 and 06.0)
Rscript --vanilla ~/bin/fig5_01.8_merge_methylation-meta-transcriptome.R
# get KOs over DMRs
Rscript --vanilla ~/bin/fig5_01.9_get-KOs-over-DMRs.R
# get the master table
Rscript --vanilla ~/bin/fig5_01.10_get-master-table.R
#####
###
# CHECK 04.4_get-DMR-positions.R I don't know what is this
# the manahttan-plots for candidate loci
