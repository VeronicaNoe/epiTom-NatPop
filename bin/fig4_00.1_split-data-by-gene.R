#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(data.table)
  library(tidyr)
  library(parallel)
})

setwd("/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_get-epialleles-over-genes/ba_DMRs-over-annotation")
outResults<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/results/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_get-epialleles-over-genes/bb_split-data-by-gene/"

geneTE<-data.table::fread("/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/gene-TE.ID", sep = '\t', data.table = T, 
                          fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
geneTE$V5<-gsub('Name=','',geneTE$V5)
geneTE<-unlist(strsplit(geneTE$V5, split = ","))

input<-list.files(pattern = ".tmp", full.names = FALSE)
#input<-grep('DMR', input, value = T, invert = T)
epialleleState <- data.table(geneName = character(),geneType = character(), strand = character(), 
                             start = integer(), annotation = character(),
                             accessions = character(), accEpiallele = character())
infoCol<-c('geneName','geneType','strand','dmr','start','DMRoverTEs','annotation')
for (i in input){
  anno<-unlist(strsplit(i, '.', fixed = T))[1]
  dmr<-unlist(strsplit(i, '.', fixed = T))[2]
  TEs<-ifelse(unlist(strsplit(i, '.', fixed = T))[4]!="tmp",
               unlist(strsplit(i, '.', fixed = T))[4], 'DMRwoTEs')
  tmp<- data.table::fread(i, sep = '\t', data.table = T, fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
  tmp$geneType <- ifelse(tmp$geneName %in% geneTE, "gene-TE", "gene")
  tmp$annotation<-anno
  tmp$dmr<-dmr
  tmp$DMRoverTEs<-TEs
  sampleCol<-grep('_leaf', colnames(tmp), value = T)
  toKeep<-c(infoCol, sampleCol)
  tmp<-tmp[,..toKeep]
  # 
  tmp<-reshape2::melt(tmp, id=infoCol,
                      variable.name="accessions", value.name = "accEpiallele")
  tmp<-as.data.table(tmp)
  epialleleState <- rbindlist(list(epialleleState, tmp), fill = TRUE)
  gc()
}
# 
# data.table::fwrite(epialleleState, file=paste0(outResults,"00.1.allSamples_epiAlleles.tsv"), quote=F,row.names=F,col.names = T,sep="\t")
#save 02
genes<-unique(epialleleState$geneName)
# 30547
#data.table::fwrite(as.data.table(genes), file=paste0(outDir,"00.2.geneList2target"), quote=F,
#                  row.names=F,col.names = F,sep="\t")
allGenes<-data.table::fread("/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/allGenes.bed", sep = '\t', data.table = T, 
                          fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
allGenes<-allGenes %>%
  filter(V1!="SL2.50ch00")
#save 03
genesWOdmrs<-setdiff(allGenes$V5, epialleleState$geneName)
data.table::fwrite(as.data.table(genesWOdmrs), file=paste0(outResults,"00.2.genes_wo_DMRs.list"), quote=F,row.names=F,col.names = T,sep="\t")

length(unique(epialleleState$geneName))
gene_split <- split(epialleleState, by = "geneName")
# Function to write each subset to a file
write_gene_file <- function(gene_data, gene_name) {
  fwrite(gene_data, file = paste0(outDir, gene_name, "_allSamples_epialleles.tsv"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
# Get the list of unique gene names
gene_names <- names(gene_split)
# Apply the function to each gene subset in parallel
mclapply(seq_along(gene_names), function(i) {
  write_gene_file(gene_split[[i]], gene_names[i])
}, mc.cores = detectCores())  
