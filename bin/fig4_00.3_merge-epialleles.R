suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(data.table)
  library(tidyr)
  library(parallel)
})

setwd("10_data-analysis/Fig4/aa_get-epialleles-over-genes/bc_gene-epialleles")
outDir<-"10_data-analysis/Fig4/aa_get-epialleles-over-genes/bd_merged-gene-epialleles/"
#### merge epialleles
# List all files with the specified pattern
inputGenes <- list.files(pattern = ".geneEpialleleState", full.names = FALSE)
inputSamples <- list.files(pattern = ".sampleEpialleleState", full.names = FALSE)
inputDMRs <- list.files(pattern = ".dmrEpialleleState", full.names = FALSE)
# Function to read a single file
read_file <- function(file) {
  data.table::fread(file, sep = '\t', data.table = TRUE, 
                    fill = TRUE, check.names = FALSE, na.strings = c("NA"), nThread = 10)
}
# Add geneLength
geneLength<-data.table::fread("05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/allGenes.bed",sep = '\t', data.table = TRUE, header = F,
                              fill = TRUE, check.names = FALSE, na.strings = c("NA"), nThread = 10)
geneLength$length<-geneLength$V3-geneLength$V2
geneLength<-geneLength[,c(5,6)]
colnames(geneLength)<-c("geneName","length")
# Detect the number of cores available
num_cores <- detectCores()
# gene epialleles
geneEpi <- mclapply(inputGenes, read_file, mc.cores = num_cores)
geneEpi <- rbindlist(geneEpi, fill = TRUE)
#30547
data.table::fwrite(geneEpi, file=paste0(outDir,"geneEpialleleState.tsv"),
                   quote=F,row.names=F,col.names = T,sep="\t")
# sample epialleles
sampleEpi <- mclapply(inputSamples, read_file, mc.cores = num_cores)
sampleEpi <- rbindlist(sampleEpi, fill = TRUE)
#30547
sampleEpi <- sampleEpi %>%
  left_join(geneLength %>% select(geneName,length), by = "geneName")
data.table::fwrite(sampleEpi, file=paste0(outDir,"samplesEpialleleState.tsv"),
                   quote=F,row.names=F,col.names = T,sep="\t")

dmrEpi <- mclapply(inputDMRs, read_file, mc.cores = num_cores)
dmrEpi <- rbindlist(dmrEpi, fill = TRUE)
data.table::fwrite(dmrEpi, file=paste0(outDir,"dmrEpialleleState.tsv"),
                   quote=F,row.names=F,col.names = T,sep="\t")
