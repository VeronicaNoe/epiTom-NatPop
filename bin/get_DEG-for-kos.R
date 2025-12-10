suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
})
setwd("/mnt/disk2/vibanez/09_KO-processing/09.1_rnaseq/ac_counts")
############### load expression counts
koList<-c("ddm1",'met1','cmt3','kyp')
input<-list.files(pattern = "_edited.tsv", full.names = FALSE)

for (k in koList){
  if(k=="kyp"){
    ko<-grep(k, input, value = T)
    wt<-grep("wt-cmt3", input, value = T)
    ko<-c(ko,wt)
  }else{
    ko<-grep(k, input, value = T) 
  }
  colNames<-c()
  rnaCounts <- c()
  for(i in ko){
    sName<-gsub('.counts_edited.tsv','',i)
    sName<-gsub('_RNAseq','',sName)
    sName<-gsub('-se','',sName)
    sName<-gsub('_mRNAseq','',sName)
    sName<-gsub('rep','',sName)
    print(sName)
    df<-data.table::fread(i, sep = '\t',
                            data.table = FALSE, fill = TRUE, check.names=FALSE,
                            na.string=c("NA"), nThread = 10)
    head(df)
    colnames(df)[2]<-sName
    toKeep<-colnames(df)[2]
    rnaCounts<-cbind(rnaCounts,df[,toKeep])
    colNames<-append(colNames,toKeep)
  }
  head(rnaCounts)
  ################ get gene names
  rnaCounts<-cbind.data.frame(rnaCounts, df[,1])
  ####### re name some samples
  colnames(rnaCounts)<-c(colNames,'ID')
  rownames(rnaCounts)<-rnaCounts$ID
  refSample<-grep('wt', colNames, value = T)
  allSample<-grep('wt', colNames, value = T, invert = T)
  rnaCounts<-rnaCounts[,c(refSample, allSample)]
  head(rnaCounts)
  myCondition<-data.frame (SampleID  = c(refSample, allSample),
                           Condition = c(rep('Control', times=length(refSample)),
                                         rep('Treat', times=length(allSample))))
  # # do the deseq transformation 
  #unadjusted
  dds <- DESeqDataSetFromMatrix(countData = rnaCounts, 
                                colData = myCondition, design = ~Condition)
  dds <- DESeq(dds)
  res<-as.data.frame(results(dds))
  res$Geneid<-rownames(res)
  head(res)
  df<-as.data.frame(counts(dds, normalized=TRUE))
  df$Geneid<-rownames(dds)
  head(df)
  mergedDf <- df %>%
    left_join(res %>% select(Geneid, log2FoldChange,padj), by = "Geneid")
  head(mergedDf)
  #save
  data.table::fwrite(mergedDf, 
                     file= paste0(k,'_counts_normalized_padj.tsv'), quote=F,
                     row.names=T,col.names = T,sep="\t")
}
