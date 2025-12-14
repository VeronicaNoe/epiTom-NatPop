#!/usr/bin/Rscript
require("data.table")
library(dplyr) 
inDir<-"05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation/"
outDir<-"06_get-meth-vcf/aa_output/"
commonDir<-"/mnt/disk2/vibanez/06_get-meth-vcf/common-accessions/data/"
setwd("06_get-meth-vcf")

infileName<-data.table::fread(paste0(inDir,"/00_colNames.tsv"), header = FALSE,sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 

args <- commandArgs(trailingOnly = TRUE)
#args<-"ch12_CG-DMR"
infoCols<-c('chr', 'start', 'end')
cName<-c(infoCols, infileName$V1)
df<- data.table::fread(paste0(inDir,args[1],'.merged.methylation.bed'), 
                       sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
colnames(df)<-cName
nameToSave<-paste0(args[1],'.tmp')
toKeep<-grep('_leaf', infileName$V1, value=TRUE)
cat(paste0("============= Getting epialleles in ", args[1], "\n"))

toSave_list <- lapply(1:nrow(df), function(r) {
  tmpDF <- df[r, toKeep, drop = FALSE]
  tmpCol <- df[r, infoCols, drop = FALSE]
  if (sum(!is.na(tmpDF)) < 2) return(NULL)
  
  densityDF <- density(unlist(tmpDF), na.rm = TRUE) #     # from https://stackoverflow.com/questions/25264461/find-local-minimum-in-bimodal-distribution-with-r
  threshold <- densityDF$x[which.min(abs(diff(densityDF$y)))]
  
  tmpDF[tmpDF > threshold] <- '1/1'
  tmpDF[!is.na(tmpDF) & tmpDF < threshold] <- '0/0'
  tmpDF[is.na(tmpDF)] <- './.'
  cbind(tmpCol, tmpDF)
})
toSave <- data.table::rbindlist(toSave_list, fill = TRUE)

sName<-grep('_leaf', colnames(toSave), value=TRUE)
#sName<-grep('-china', sName, value=TRUE, invert = TRUE)
sample2Remove<-c("LA0534_leaf","Ohio-8245_leaf")
sName<-setdiff(sName,sample2Remove)

infoCol<-toSave[,c("chr","start")]
samplesCol<-toSave[,..sName]
#head(samplesCol)

cat(paste0("============= Formating for plink sample ", args[1], "\n"))
REF<-data.frame(samplesCol[,"TS-253_leaf"])
ALT<-REF
colnames(ALT)<-"ALT"

ALT<-gsub("0", "C", ALT$ALT) # just to get the alternative
ALT<-gsub("1", "T", ALT)
ALT<-as.data.frame(ALT)

colnames(REF)<-"REF"
REF<-gsub("0", "T", REF$REF)	# 0 is not methylated
REF<-gsub("1", "C", REF)	# 1 is methylated
REF<-as.data.frame(REF)

ID<-paste0(toSave$chr,":",toSave$start)
QUAL<-rep("30", times=nrow(samplesCol))
FILTER<-rep(".", times=nrow(samplesCol))
INFO<-rep("NA", times=nrow(samplesCol))
FORMAT<-rep("GT", times=nrow(samplesCol))

cat(paste0("============= Saving sample ", args[1], "\n"))

out.vcf<-cbind.data.frame(infoCol, ID, REF, ALT, QUAL, FILTER, INFO,FORMAT, samplesCol)
colnames(out.vcf)[c(1,2)]<-c("#CHROM","POS")
#head(out.vcf)
# write.table(out.vcf, paste0(commonDir,nameToSave),quote=F,row.names=F,sep="\t")
data.table::fwrite(out.vcf, file=paste0(commonDir,nameToSave),quote=F,row.names=F,col.names = T,sep="\t", nThread=4)
#
# save data to 
sName<-grep('-china', sName, value=TRUE, invert = TRUE)
tmp<-samplesCol[,..sName]
out.vcf<-cbind.data.frame(infoCol, ID, REF, ALT, QUAL, FILTER, INFO,FORMAT, tmp)
colnames(out.vcf)[c(1,2)]<-c("#CHROM","POS")
#write.table(out.vcf, paste0(outDir,nameToSave),quote=F,row.names=F,sep="\t")
data.table::fwrite(out.vcf, file=paste0(outDir,nameToSave),quote=F,row.names=F,col.names = T,sep="\t", nThread=4)

### split by groups and annotations
anno<-c("gene-TE","gene","intergenic","TE-wo-gene","promotor")
group<-c("all","1-WILD","2-PIM","3-SLC","4-SLL")
dmr<-unlist(strsplit(args[1],'_',fixed = T))[2]
infoVcf<-grep('_leaf',colnames(out.vcf), invert = T, value = T)
for(a in anno){
  for(g in group){
  anno_dir<-"05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection/"
  anno_lis<-data.table::fread(paste0(anno_dir,a,".",dmr,'.methylation'), 
                           sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
  colnames(anno_lis)<-cName
  anno_lis$ID<-paste0(anno_lis$chr,":",anno_lis$start)
  anno.vcf <- out.vcf %>%
    filter(ID %in% anno_lis$ID)
  colnames(anno.vcf)[c(1,2)]<-c("#CHROM","POS")
  if(g=="all"){
    data.table::fwrite(out.vcf, file=paste0(outDir,args[1],"_",a,".tmp"),quote=F,row.names=F,col.names = T,sep="\t", nThread=4)  
  }else{
    group_list<-data.table::fread(paste0(outDir,g,'.group'), header = F,
                                  sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
    anno.vcf <- anno.vcf %>%
      select(infoVcf,group_list$V1)
    data.table::fwrite(out.vcf, file=paste0(outDir,args[1],"_",a,"_",g,".tmp"),quote=F,row.names=F,col.names = T,sep="\t", nThread=4)  
  }
  }
}
