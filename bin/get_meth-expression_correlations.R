#!/usr/bin/Rscript
suppressPackageStartupMessages({
  library(tidyr)
})
args <- commandArgs(trailingOnly = TRUE)
#set paths and chr

setwd("10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/bc_format-data")
outDir<-"10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/bd_get-DMR-meth-gene-expression-correlation/tmp/"

print("###### Loading meth data")
chr<-args[1]
dmr<-args[2]

methLevels<-data.table::fread(paste0(chr,'_',dmr,'_methylation-levels.tsv'), sep = '\t',
                              data.table = FALSE, fill = TRUE, check.names=FALSE,
                              na.string=c("NA"), nThread = 10)
toKeep<-args[3]

methLevels$geneName <- sub("\\.\\d+\\.\\d+$", "", methLevels$geneName)

methLevels<-subset(methLevels, toFilter==toKeep)
gene<-unique(methLevels$geneName)

print("###### Loading gene data")

geneExpression<-data.table::fread('07_rnaseq-processing/07.3_counts/leaf_transcriptome-counts_normalized_batch-adjusted.tsv', sep = '\t',
                                  data.table = FALSE, fill = TRUE, check.names=FALSE,
                                  na.string=c("NA"), nThread = 10)

print("###### Calculating correlations")
for(g in gene){
  tmpMethLevels<-subset(methLevels, geneName==g)
  tmpGeneExpression<-subset(geneExpression, Geneid==g)
  tmpGeneExpression<-tmpGeneExpression %>%
    pivot_longer(
      cols = -Geneid,  # Exclude the Geneid column
      names_to = "Sample",  # Name for the column holding previous column names
      values_to = "Expression"  # Name for the column holding values
    )
  tmpGeneExpression$expressionLog2<-log2(tmpGeneExpression$Expression+1)
  #head(geneExpression)
  tmp <- dplyr::left_join(tmpMethLevels, tmpGeneExpression, by = "Sample")
  # Perform Pearson correlation test
	result <- tryCatch({
		pearsonExpMeth <- cor.test(tmp$expressionLog2, tmp$methylationLevels, na.rm = TRUE)
		pValpearsonExpMeth <- pearsonExpMeth$p.value
		pearsonExpMeth <- pearsonExpMeth$estimate
		list(pValue = pValpearsonExpMeth, estimate = pearsonExpMeth)
		}, warning = function(w) {
		pValpearsonExpMeth <- "NA"
		pearsonExpMeth <- 0
		list(pValue = pValpearsonExpMeth, estimate = pearsonExpMeth)
		}, error = function(e) {
		pValpearsonExpMeth <- "NA"
		pearsonExpMeth <- 0
		list(pValue = pValpearsonExpMeth, estimate = pearsonExpMeth)
	})
	pValpearsonExpMeth <- result$pValue
	pearsonExpMeth <- result$estimate
  # Perform Spearman correlation test
	result <- tryCatch({
		spearmanExpMeth <- cor.test(tmp$expressionLog2, tmp$methylationLevels, method = 'spearman', na.rm = TRUE)
		pValSpearmanExpMeth <- spearmanExpMeth$p.value
		spearmanExpMeth <- spearmanExpMeth$estimate
		list(pValue = pValSpearmanExpMeth, estimate = spearmanExpMeth)
		}, warning = function(w) {
		pValSpearmanExpMeth <- "NA"
                spearmanExpMeth <- 0
                list(pValue = pValSpearmanExpMeth, estimate = spearmanExpMeth)
		}, error = function(e) {
		pValSpearmanExpMeth <- "NA"
		spearmanExpMeth <- 0
                list(pValue = pValSpearmanExpMeth, estimate = spearmanExpMeth)
	})
	pValSpearmanExpMeth <- result$pValue
	spearmanExpMeth <- result$estimate
#check
	if(length(unique(tmp$genePart))>1){
	  genePart<-paste0(unique(tmp$genePart), collapse = ",")
	}else{
	  genePart<-unique(tmp$genePart)
	}

  corOut <- data.frame(matrix(nrow = 1, ncol = 10))
  corOut[1,1]<-unique(tmp$chr)
  corOut[1,2]<-unique(tmp$start)
  corOut[1,3]<-unique(tmp$toFilter)
  corOut[1,4]<-unique(tmp$Geneid)
  corOut[1,5]<-unique(tmp$strand)
  corOut[1,6]<-genePart
  corOut[1,7]<-pearsonExpMeth
  corOut[1,8]<-pValpearsonExpMeth
  corOut[1,9]<-spearmanExpMeth
  corOut[1,10]<-pValSpearmanExpMeth

  colnames(corOut)<-c('chr','start','DMR','geneID','strand','genePart','pearsonExpMeth','pValpearsonExpMeth','spearmanExpMeth','pValSpearmanExpMeth')

print("################ SAVING correlation")
  data.table::fwrite(corOut,
                   file= paste0(outDir,toKeep,"_",g,"_expression.correlation"), quote=F,
                   row.names=F,col.names = T,sep="\t")
}

