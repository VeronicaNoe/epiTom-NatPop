#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
path2results<-"/mnt/disk2/vibanez/06_DMR/tmp/"
#path2acc<-"/mnt/disk2/vibanez/05_methylkit/results/DMR/chr01/"
#toTest<-"P12"
#toUnique<-"_DMR
accFactors<-unlist(strsplit(args[1],"_"))
accName<-accFactors[2]
chr2process=accFactors[1]
cat(accName)
meth<-paste0(accName,'_methylated')
unmeth<-paste0(accName,'_unmethylated')
cnames<-c(meth, unmeth)
  
toKeep<-c('chr', 'start', 'strand', cnames, 'context')  
ctxt<-c("CG", "CHG", "CHH")
DMR<- data.table::fread(args[1], sep = '\t', data.table = FALSE, fill = TRUE,header = T, na.string=c("NA"), nThread = 20) 

for (c in ctxt){
  subDMR<-subset(DMR, context==c)
  subDMR<-subDMR[,toKeep]
  data.table::fwrite(subDMR, file=paste0(path2results, chr2process,"_",c,"_",accName,"_DMR.bed"), quote=F,row.names=F,col.names = F,sep="\t")
}
