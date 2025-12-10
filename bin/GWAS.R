#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
# get the folder passed from the shell
suppressPackageStartupMessages({
  require("data.table")
  library(ggplot2)
  #library(qqman)
})
args <- commandArgs(trailingOnly = TRUE)
setwd("/home/vibanez/results")
print("##########  Loading sample files  ###########")
samples<-list.files(pattern = ".gz", full.names = FALSE)
FROM<-args[1]
TO<-args[2]

if(TO>length(samples)){
  TO<-length(samples)}

mQTL<-c()
result<-c()
print("Merging data")
for(i in FROM:TO){
	if(length(readLines(samples[i], warn = FALSE))==0){
      next }
	df<- data.table::fread(samples[i], sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
	chr<-unlist(strsplit(samples[i],"_", fixed = TRUE))[1]
	dmr<-unlist(strsplit(samples[i],"_", fixed = TRUE))[2]
	pstn<-unlist(strsplit(samples[i],"_", fixed = TRUE))[3]
	pstn<-gsub(".filtered.ps.gz", "", pstn)
	df<-subset(df, V4<5e-8)
	if(length(readLines(paste0(chr,"_",dmr,"_",pstn,".reml"), warn = FALSE))==0){
		next }
	reml<-t(data.table::fread(paste0(chr,"_",dmr,"_",pstn,".reml"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20))
	if(dim(df)[1]==0){
		rTemp<-cbind.data.frame(chr,dmr,pstn,"0", reml)
		colnames(rTemp)[4]<-'nrow(df)'
		result<-rbind.data.frame(result, rTemp)
		next
	} else {
		Temp<-cbind.data.frame(chr,dmr,pstn,nrow(df), reml)
		result<-rbind.data.frame(result, rTemp)
		snpInfo<-do.call(rbind, strsplit(as.character(df$V1), ":", fixed = TRUE))
		DMR<-rep(dmr, times=nrow(df))
		dmrPstn<-rep(pstn, times=nrow(df))
		dmrChr<-rep(gsub("ch","",chr), times=nrow(df))
		mQTLtmp<-cbind.data.frame(snpInfo,df[,2:4],DMR,dmrChr, dmrPstn)
		mQTL<-rbind.data.frame(mQTL,mQTLtmp)
	}
}
colnames(result)<-c('chr', 'dmr', 'position', 'mQTL',"logVariance",
                    "log", "delta", "sigmaGen", "sigmaError", "h2")
# #camus
outDir<-c("/mnt/disk2/vibanez/GWAS/analysis-Results/")
print("Saving data")
data.table::fwrite(result, file=paste0(outDir, as.character(FROM), "_summaryGWAS.tsv"), quote=F,row.names=F,col.names = T,sep="\t")
if(length(mQTL)>0){
	colnames(mQTL)<-c('chrSNP','startSNP','Beta','SEbeta','pVal','DMR','chrDMR','startDMR')
	data.table::fwrite(mQTL, file=paste0(outDir, as.character(FROM), "_mQTL.tsv"), quote=F,row.names=F,col.names = T,sep="\t")
}
# #mcclintock
#data.table::fwrite(result, file="/mnt/data6/vibanez/SNPs/vcfiles/ah_GWAS/summaryGWAS.tsv", quote=F,row.names=F,col.names = T,sep="\t")
#data.table::fwrite(mQTL, file="/mnt/data6/vibanez/SNPs/vcfiles/ah_GWAS/mQTL.tsv", quote=F,row.names=F,col.names = T,sep="\t")
