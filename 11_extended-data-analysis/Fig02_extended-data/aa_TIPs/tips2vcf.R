setwd("/home/IPS2/vibanez/Desktop/Q-lab/2023/TIPs")
# load all sig and nonsig
tips<-data.table::fread('TIPs.csv', sep = '\t', data.table = FALSE, 
                              fill = TRUE, nThread = 20)
tips[1:6,1:6]
tips$TE<-gsub('_','-', tips$TE)
tmpCol <- do.call(rbind, strsplit(tips$TE, "-"))
tmp <- data.frame(tmpCol, stringsAsFactors = FALSE)
head(tmp)

sampleNames<-data.table::fread('sample2keep', sep = '\t', data.table = FALSE, 
                        fill = TRUE,header = F, nThread = 20)
head(sampleNames)
rownames(sampleNames)<-sampleNames$V1
toKeep<-intersect(sampleNames$V1,colnames(tips))
sampleNames<-sampleNames[toKeep,]
infoCol<-tips[,c("chr","start")]
infoCol$chr<-gsub('SL2.50ch0','',infoCol$chr)
samplesCol<-tips[,toKeep]
colnames(samplesCol)<-sampleNames$V2
samplesCol[samplesCol ==1]<-"1|1"
samplesCol[samplesCol ==0 ]<-"0|0"
samplesCol[samplesCol =="NA" ]<-".|."
head(samplesCol)

REF<-data.frame(samplesCol[,"TS-253"])
ALT<-REF
colnames(ALT)<-"ALT"

ALT<-gsub("1", "T", ALT$ALT)
ALT<-gsub("0", "C", ALT)
ALT<-as.data.frame(ALT)

colnames(REF)<-"REF"
REF<-gsub("0", "T", REF$REF)
REF<-gsub("1", "C", REF)
REF<-as.data.frame(REF)

ID<-paste0(tips$chr,":",tips$start,":",tmp$X1)
QUAL<-rep("30", times=nrow(samplesCol))
FILTER<-rep(".", times=nrow(samplesCol))
INFO<-rep("NA", times=nrow(samplesCol))
FORMAT<-rep("GT", times=nrow(samplesCol))
tipsvcf<-cbind.data.frame(infoCol, ID, REF, ALT, QUAL, FILTER, INFO,FORMAT, samplesCol)
colnames(tipsvcf)[c(1,2)]<-c("#CHROM","POS")
head(tipsvcf)

#### add SNPs
snps<- data.table::fread("SNPs-TIP.LD.edited.vcf",
                         sep = '\t', data.table = FALSE, fill = TRUE, na.string="NA", nThread = 20) 
snpNames<- data.table::fread("sample4tips",sep = '\t', data.table = FALSE,
                             header = F,
                             fill = TRUE, na.string="NA", nThread = 20) 
infoCol<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
colnames(snps)<-c(infoCol, snpNames$V1)
head(snps)
#
tipSamples<-setdiff(colnames(tipsvcf), infoCol)
snpSamples<-setdiff(colnames(snps), infoCol)
sameSamples<-intersect(tipSamples,snpSamples)
tipsvcf<-tipsvcf[,c(infoCol,sameSamples)]
data.table::fwrite(tipsvcf, file= "TIPs.vcf", quote=F,
                   row.names=F,col.names = T,sep="\t")

snps<-snps[,c(infoCol,sameSamples)]
# merge
snps[1:6,1:10]
tipsvcf[1:6,1:10]

out<-rbind.data.frame(tipsvcf,snps)
out<-out[order(out$`#CHROM`),]

data.table::fwrite(out, file= "TIPs-SNPs.vcf", quote=F,
                   row.names=F,col.names = T,sep="\t")
