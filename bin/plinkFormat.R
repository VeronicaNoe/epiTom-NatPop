#!/users/bioinfo/vibanez/miniconda3/envs/methylKit/bin/Rscript
require("data.table")
args <- commandArgs(trailingOnly = TRUE)
setwd("/mnt/disk2/vibanez/02_methylkit/af_mVCF/aa_data")
output<-"/mnt/disk2/vibanez/02_methylkit/af_mVCF/ab_mvcf-files/"

input<-list.files(pattern = ".meth4plink", full.names = FALSE)
dfNames<- data.table::fread("/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/bb_output/00_colNames.tsv", header = FALSE,sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
infileName<-gsub("_Biseq_CG_chrNonRef.methylation.bed", "",dfNames[,1])

anno<-args[2]
chr<-args[3]
DMR<-args[4]
threshold<-args[5]

#df<- data.table::fread("ch01_gene-wo-TE.C-DMR.meth4plink", sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
df<- data.table::fread(args[1], sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
colnames(df)<-c("chr", "start","end",infileName)
#toKeep<-intersect(sampleTab$sample, colnames(df)) # remove fruit
#df<-df[,c("chr","start","end", toKeep)]
infoCol<-df[,c("chr","start")]
#samplesCol<-df[,toKeep]
samplesCol<-df[,infileName]

samplesCol[samplesCol >=threshold]<-"1|1"
samplesCol[samplesCol !=".|." & samplesCol < threshold ]<-"0|0"

cat(paste0("============= Processing sample ", args[1], "\n"))
REF<-data.frame(samplesCol[,"TS-253_leaf"])
colnames(REF)<-"TS-253_leaf"


ALT<-REF
ALT<-gsub("1", "T", ALT$`TS-253_leaf`)
ALT<-gsub("0", "C", ALT)
ALT<-as.data.frame(ALT)
colnames(ALT)<-"ALT"

REF<-gsub("0", "T", REF$`TS-253_leaf`)
REF<-gsub("1", "C", REF)
REF<-as.data.frame(REF)
colnames(REF)<-"REF"

ID<-paste0(df$chr,df$start)
QUAL<-rep("30", times=nrow(samplesCol))
FILTER<-rep(".", times=nrow(samplesCol))
INFO<-rep("NA", times=nrow(samplesCol))
FORMAT<-rep("GT", times=nrow(samplesCol))

cat(paste0("============= Saving sample ", args[1], "\n"))
out.vcf<-cbind.data.frame(infoCol, ID, REF, ALT, QUAL, FILTER, INFO,FORMAT, samplesCol)
colnames(out.vcf)[c(1,2)]<-c("#CHROM","POS")
#out<-rbind("##fileformat=VCFv4.3", out.vcf)
write.table(out.vcf, paste0(output,chr,"_",anno,"_",DMR,".tmp"),quote=F,row.names=F,sep="\t")

