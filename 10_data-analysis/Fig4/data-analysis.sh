#!/bin/bash
## get the epiallele state of samples along with its gene information
Rscript --vanilla ~/bin/fig4_01.1_merge_geneInfo-sampleInfo-TEdistance.R

## get the epiallele frequency over genes and tomato groups
Rscript --vanilla ~/bin/fig4_01.2_epiallele-frequency_over-genes.R

## get the epiallele effect of having TEs sourrounding the gene
Rscript --vanilla ~/bin/fig4_01.3_TEs-effect_over-gene-epialleles.R

## get gene expression per each epiallele gene
Rscript --vanilla ~/bin/fig4_01.4_gene-expression_epiallele.R

## get genes that are DE in kos and also when they have different meth state
Rscript --vanilla ~/bin/fig4_01.5_merge_epialleles-DEG_with_KOs-DEG.R

### get geneList coordinates
ANNO_DIR="05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
RESULT_DIR="10_data-analysis/Fig4/results"
grep -Ff ${RESULT_DIR}/gene_list.tsv ${ANNO_DIR}/allGenes.bed > ${RESULT_DIR}/geneList_coordinates.bed

