#!/usr/bin/Rscript
suppressPackageStartupMessages({
	library(dplyr)
	library(plyr)
	library(data.table)
	library(tidyr)
})

setwd("10_data-analysis/Fig4/aa_get-epialleles-over-genes/bb_split-data-by-gene")
outDir<-"10_data-analysis/Fig4/aa_get-epialleles-over-genes/bc_gene-epialleles/"
args <- commandArgs(trailingOnly = TRUE)

epialleleState<- data.table::fread(paste0(args[1],"_allSamples_epialleles.tsv"), sep = '\t', 
                                   data.table = T, fill = TRUE, check.names=FALSE,
                                   na.string=c("NA"), nThread = 20)
epialleleState$annotation<-gsub("-TE",'',epialleleState$annotation)
define_sampleEpiallele <- function(sumSampleEpi, sumSampleDMR) {
  # Handle cases where sumSampleEpi is NA or contains invalid values
  if (is.na(sumSampleEpi) || sumSampleEpi %in% c("NA", "0", "0,NA", "NA,0")) {
    return("UM")
  } else if (sumSampleEpi %in% c("1", "0,1", "1,0", "1,NA", "NA,1")) {
    if (sumSampleDMR == "CG-DMR") {
      return("gbM")
    } else if (sumSampleDMR == "C-DMR") {
      return("teM")
    } else if (sumSampleDMR == "C-DMR,CG-DMR") {
      if (sumSampleEpi %in% c("1,0", "1,NA", "1")) {
        return("teM")
      } else if (sumSampleEpi %in% c("0,1", "NA,1")) {
        return("gbM")
      }
    }
  }
  return(NA)
}

cat("### get dmr Epiallele\n")
dmrEpi <- epialleleState %>%
  filter(geneName==args[1]) %>%
  group_by(geneName, geneType, strand, start, annotation, accessions, DMRoverTEs) %>%
  dplyr::summarise(
    sumSampleEpi = paste(unique(accEpiallele), collapse = ","),
    sumSampleDMR = paste(unique(dmr), collapse = ","),
    dmrEpiallele = define_sampleEpiallele(sumSampleEpi, sumSampleDMR),
    .groups = 'drop'
  ) %>%
  as.data.table()
gc()
data.table::fwrite(dmrEpi, file=paste0(outDir,args[1],".dmrEpialleleState"), quote=F,row.names=F,col.names = T,sep="\t")
cat("### get sample Epiallele\n")
sampleEpi <- dmrEpi %>%
    group_by(geneName, geneType, strand, annotation, accessions) %>%
    dplyr::summarise(
      sumSampleEpi = paste(unique(dmrEpiallele), collapse = ","),
      UMn= sum(dmrEpiallele=="UM"),
      teMn= sum(dmrEpiallele=="teM"),
      gbMn= sum(dmrEpiallele=="gbM"),
      nDMR= sum(UMn,teMn,gbMn),
      sampleEpiallele = ifelse(teMn !=0 & teMn > UMn & teMn > gbMn, "teM",
                               ifelse(gbMn !=0 & gbMn > UMn & gbMn >= teMn, "gbM","UM")),
      teMnoverTEs = sum(ifelse(DMRoverTEs == "DMRoverTEs" & sampleEpiallele == "teM", 1, 0)),
      .groups = 'drop')

cat("### get gene Epiallele\n")
geneEpi<-sampleEpi %>%
  group_by(geneName, strand, geneType, annotation) %>%
  dplyr::summarise(
    sumSampleEpi = paste(unique(sampleEpiallele), collapse = ","),
    nDMR = round((sum(gbMn)+sum(teMn)+sum(UMn))),
    geneEpiallele = ifelse(grepl("teM", sumSampleEpi) & grepl("gbM", sumSampleEpi), "PE", 
			ifelse(grepl("gbM", sumSampleEpi), "gbM",
				ifelse(grepl("teM", sumSampleEpi), "teM", 
                                         "UM"))),
    gbMn=sum(gbMn),
    teMn=sum(teMn),
    UMn=sum(UMn),
    .groups = 'drop')
data.table::fwrite(geneEpi, file=paste0(outDir,args[1],".geneEpialleleState"), quote=F,row.names=F,col.names = T,sep="\t")
gene<- geneEpi %>%
  group_by(geneName, geneType) %>%
  dplyr::summarise(
    sumGeneEpiallele = paste(unique(geneEpiallele), collapse = ","),
    sumAnnotation = paste(unique(annotation), collapse = ","),
    var= ifelse(sumGeneEpiallele  %in% c("gbM","teM","PE") & sumAnnotation=="exon,intron", "variable",
                "nonVar"))

sampleEpi <- sampleEpi %>%
  left_join(geneEpi %>% select(geneName,annotation ,geneEpiallele), by = c("geneName", "annotation"))
sampleEpi <- sampleEpi %>%
  left_join(gene %>% select(geneName,geneType, var), by = c("geneName", "geneType"))
data.table::fwrite(sampleEpi, file=paste0(outDir,args[1],".sampleEpialleleState"), quote=F, row.names=F,col.names = T,sep="\t")
