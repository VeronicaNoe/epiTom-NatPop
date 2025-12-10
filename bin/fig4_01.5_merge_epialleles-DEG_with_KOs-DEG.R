#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggpubr)
  library(plyr)
  library(ggplot2)
  library(data.table)
  library("RColorBrewer")
  library("gplots")
  library(viridis)
  library(tidyr)
})
#in Camus:/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/toMergeFiles
# setwd("~/Desktop/Q-lab/IBANEZ_etal_2024/Fig4/data")
# outDir<-"~/Desktop/Q-lab/IBANEZ_etal_2024/Fig4/results/"
# outPlot<-"~/Desktop/Q-lab/IBANEZ_etal_2024/Fig4/plots/"
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_get-epialleles-over-genes/bd_merged-gene-epialleles")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/plots/"


########################################################################################
# merge all the information to get:
# geneName | geneType | geneEpiallele | annotationMerged | strand | absDistanceTE | 
# TEposition | accessions |  accGroup | sampleEpiallele
########################################################################################
## all genes with consistent epialleles: 
# outDir,'00_epiallele_sample.tsv'
## genes with variable epialleles: 
# outDir,'00_varaible-epiallele-overGenes-w-sample.tsv'
#### in camus
# add group info
sampleTab<-data.table::fread("00_sampleTab.tsv", header = TRUE,sep = '\t',
                             data.table = T, fill = TRUE, na.string=c("NA"), nThread = 20)
sampleTab<-subset(sampleTab, Organ=='leaf')
# get epialle states for samples
epialleleState<- data.table::fread("samplesEpialleleState.tsv", sep = '\t', 
                                   data.table = T, fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
epialleleState$accessions<-gsub('_leaf','',epialleleState$accessions)
unique(epialleleState$sampleEpiallele)

# get TE distance to genes
geneTEdistance<- data.table::fread("/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/TE_closest_gene.tsv", sep = '\t', 
                                   data.table = T, fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
colnames(geneTEdistance)[4]<-"strand"
colnames(geneTEdistance)[5]<-"geneName"
colnames(geneTEdistance)[6]<-"teDistance"
geneTEdistance <- geneTEdistance %>%
  filter(V1 != "SL2.50ch00") %>%
  mutate(
    absDistance = abs(teDistance),  # Absolute value of distance
    tePos = case_when(
      strand == "+" & teDistance < 0 ~ "upstream",    
      strand == "+" & teDistance > 0 ~ "downstream",  
      strand == "-" & teDistance < 0 ~ "downstream",  
      strand == "-" & teDistance > 0 ~ "upstream",    
      TRUE ~ "over-gene"                            
    )
  )

geneTEdistance
geneTEdistance <- geneTEdistance %>%
  separate_rows(geneName, sep = ",")
gc()

### merge
mergedEpiallele <- epialleleState %>%
  left_join(sampleTab[, c("BiseqName", "Group")], by = c("accessions" = "BiseqName"))
mergedEpiallele <- mergedEpiallele %>%
  left_join(geneTEdistance[, c("geneName","absDistance", "tePos")], by = c("geneName" = "geneName"))
gc()
mergedEpiallele
totalGenes<-unique(mergedEpiallele$geneName)
length(totalGenes) #30547

sumES <- mergedEpiallele[, .(
  geneEpiallele = paste(unique(geneEpiallele), collapse = ","),
  annotation = paste(unique(annotation), collapse = ",")
), by = .(geneName, geneType)]

varGenes<-subset(sumES, geneEpiallele=="PE,teM" | geneEpiallele=="PE,gbM")
sort(unique(sumES$geneEpiallele))
# Filter and modify subsets for specific allele patterns
sameAllele <- c("gbM", "teM", "PE","UM")
sameAlleleDT <- sumES[geneEpiallele %in% sameAllele]
subset(sumES, geneName=="Solyc01g005170.1")

intronPE <- c("gbM,PE", "teM,PE", "UM,PE")
intronPEdt <- sumES[geneEpiallele %in% intronPE]
intronPEdt[, geneEpiallele := "PE"]
intronPEdt[, annotation := "intron"]

exonPE <- c("PE,gbM", "PE,teM","PE,UM")
exonPEdt <- sumES[geneEpiallele %in% exonPE]
exonPEdt[, geneEpiallele := "PE"]
exonPEdt[, annotation := "exon"]

gb <- c("gbM,UM", "UM,gbM","gbM,PE","PE,gbM")
gbdt <- sumES[geneEpiallele %in% gb]
gbdt[, annotation := ifelse(geneEpiallele == "gbM,UM", "exon",
                            ifelse(geneEpiallele == "UM,gbM", "intron", 
                                   ifelse(geneEpiallele == "gbM,PE", "exon", 
                                          "intron")))]
gbdt[, geneEpiallele := "gbM"]

te <- c("teM,UM", "UM,teM","teM,PE","PE,teM")
tedt <- sumES[geneEpiallele %in% te]
tedt[, annotation := ifelse(geneEpiallele == "teM,UM", "exon", 
                            ifelse(geneEpiallele == "UM,teM", "intron", 
                                   ifelse(geneEpiallele == "teM,PE", "exon", 
                                          "intron")))]
tedt[, geneEpiallele := "teM"]
#PE genes in both, intron and exon
tmpVariablePE <- c("teM,gbM", "gbM,teM")
tmpVariablePEdt <- sumES[geneEpiallele %in% tmpVariablePE]
tmpVariablePEdt[, annotation := ifelse(geneEpiallele == "teM,gbM", "exon", "intron")]
tmpVariablePEdt[, geneEpiallele := "teM"]

variablePEdt <- data.table()
variablePEdt <- rbindlist(list(variablePEdt, tmpVariablePEdt), fill = TRUE)
variableGeneList<-unique(variablePEdt$geneName)
# Combine all subsets into one data table for this iteration
byParts <- data.table()
byParts <- rbindlist(list(byParts, sameAlleleDT, intronPEdt, exonPEdt,tedt,gbdt), fill = TRUE)
byParts <- byParts %>%
  filter(!geneName %in% varGenes$geneName) 
byParts <- byParts %>%
  separate_rows(annotation, sep = ",")
geneList<-unique(byParts$geneName)
#compare
length(totalGenes)
#[1] 30547
length(varGenes$geneName)
# 852
length(geneList)
#26279
length(variableGeneList)
#[1] 3416
data.table::fwrite(byParts,paste0(outDir,'01.0_geneEpialleles_exon-intron.tsv'),
                   quote=F, row.names=F,
                   col.names = T,sep="\t")

data.table::fwrite(variablePEdt,paste0(outDir,'01.1_geneEpialleles_variable-exon-intron.tsv'),
                   quote=F, row.names=F,
                   col.names = T,sep="\t")

mergedEpiallele<- merge(byParts,mergedEpiallele,
                        by = c("geneName", "annotation"), all.x = TRUE,
                        allow.cartesian=TRUE)%>%
  select(-geneType.y,-geneEpiallele.y, -var,-sumSampleEpi )
setDT(mergedEpiallele)
# 
# # get the epiallele of the accessions
colnames(mergedEpiallele)[2]<-"annotation"
colnames(mergedEpiallele)[3]<-"geneType"
colnames(mergedEpiallele)[4]<-"geneEpiallele"
colnames(mergedEpiallele)[13]<-"geneLength"
colnames(mergedEpiallele)[14]<-"accGroup"
colnames(mergedEpiallele)[15]<-"absDistanceTE"
colnames(mergedEpiallele)[16]<-"TEposition"
mergedEpiallele
gc()
# subset(mergedEpiallele, geneName=="Solyc01g005170.1")
# #save
data.table::fwrite(mergedEpiallele,paste0(outDir,'01.2_epiallele_sample.tsv'),
                   quote=F, row.names=F,
                   col.names = T,sep="\t")
# gc()
# unique(mergedEpiallele$geneEpiallele)
