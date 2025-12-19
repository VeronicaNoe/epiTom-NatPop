setwd("05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation")
inMeth<-list.files(pattern = ".methylation", full.names = FALSE)
methLevels<-c()
for(m in inMeth){
  infoCol<-c('chr','start', 'end')
  dmr<-unlist(strsplit(m, '.', fixed=TRUE))[2]
  anno<-unlist(strsplit(m, '.', fixed=TRUE))[1]
  df<-data.table::fread(m, sep = '\t', data.table = TRUE,
                        fill = TRUE, na.string=c("NA"), nThread = 20)
  colNames<-data.table::fread(paste0(dmr,'_colNames.tsv'), sep = '\t',
                              data.table = FALSE,header = FALSE,
                              fill = TRUE, na.string=c("NA"), nThread = 20)
  colnames(df)<-c(infoCol, colNames$V1)
  DMRs<-paste0(gsub('SL2.50','',df$chr),'_',dmr,'_',df$start)
  df<-cbind.data.frame(df,anno, DMRs)
  methLevels<-rbind.data.frame(methLevels,df)
}
head(methLevels)
data.table::fwrite(methLevels, file="DMR_methylationLevels" , quote=F,
                   row.names=F,col.names = T,sep="\t", nThread = 60)
