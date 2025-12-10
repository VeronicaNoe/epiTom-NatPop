#!/usr/bin/R
require("data.table")
#library(ggplot2)
library(plyr)
#WD<-"/mnt/disk2/vibanez/GWAS"
WD<-"/mnt/disk2/vibanez/GWAS/newColab"
setwd(WD)
methDir<-paste0(WD,"/methFiles")
outDir<-paste0(WD,"/phenotype")
# get the folder passed from the shell
args <- commandArgs(trailingOnly = TRUE)
# go to within annotation folder
# load data
cat(paste0("========= Loading ", args[1], "\n"))
fileNames<-unlist(strsplit(args[1],".",fixed = TRUE))[1]
chr<-unlist(strsplit(fileNames,"_", fixed = TRUE))[1]
dmr<-unlist(strsplit(fileNames,"_", fixed = TRUE))[2]
print(paste0(chr, " ", dmr))
df<- data.table::fread(paste0(methDir,"/",args[1]), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
df$V109<-NULL
# get the fam names
sampleNames<-data.table::fread(paste0(WD,"/data/SNPs_SL2-5.LD.tfam"), 
sep = ' ', data.table = FALSE, fill = TRUE,header = FALSE, na.string=c("NA"), nThread = 20) 
methNames<-data.table::fread(paste0(methDir,"/","00_colNames.tsv"), 
sep = '\t', data.table = FALSE, fill = TRUE,header = FALSE, na.string=c("NA"), nThread = 20) 
# change name of this samples
methNames<-gsub("S-pimBGV006775","BGV006775",methNames[,1])
methNames<-gsub("LA2838A","S-lycLA2838A",methNames)
methNames<-gsub('_leaf', '', methNames)
infoCol<-c("chr","start", "end")
colnames(df)<-c(infoCol,methNames)
# change SNPs names to match with methNames
toKeep<-gsub("_1", "",sampleNames[,1])
toKeep<-gsub("S_","S-", toKeep)
# get the intersection and sort samples according to FAM files
methNames<-colnames(df)[! colnames(df) %in% infoCol]
toSort<-c(infoCol,intersect(toKeep, methNames))
out<-df[,toSort]
#head(out)
# keep SNPs names
colnames(out)<-c(infoCol, sampleNames[,1])
# format the final output
out$end<-NULL
out$chr<-NULL
rownames(out)<-out$start
out$start<-NULL
out<-t(out)
out <- data.frame(INDID = row.names(out), out, check.names = FALSE)
out <- data.frame(FAMID = row.names(out), out, check.names = FALSE)
# save it
toTarget<-c()
#toStop<-ncol(out)
#toStop<-45000
#batch<-30000
#i<-3
#while (batch<= toStop){
#save each DMR as a separate file
for (c in 3:ncol(out)){
	#print(paste0("", batch))
	sample2process<-paste0(chr,"_",dmr,"_",colnames(out)[c],".phenotype")
	toTarget<-rbind(toTarget,sample2process)
	write.table(out[,c(1,2,c)], paste0(outDir,"/",sample2process),quote=F,row.names=F,sep="\t", col.names=FALSE)
}
#save a file with all the targets
write.table(toTarget[,1], paste0(WD,"/",chr,"_",dmr,".toTarget"),quote=F,row.names=F,sep="\t", col.names=FALSE)
#	i<-batch+1
#	batch<-batch+30000
#	print(paste0("========== batch ",batch))
#}
