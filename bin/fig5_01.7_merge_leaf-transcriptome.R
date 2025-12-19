suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(DESeq2)
    library(sva)
  })
  
  setwd("07_rnaseq-processing/07.3_counts")
  outDir<-"10_data-analysis/Fig5/ab_data-analysis/results/"
  outPlot<-"10_data-analysis/Fig5/ab_data-analysis/plots/"
  
  input<-list.files(pattern = "_edited.tsv", full.names = FALSE)
  input<-grep("TS-179",input,value = T ,invert = T)
  colNames<-c()
  rnaCounts <- c()
  for (i in input){
      df<-data.table::fread(i, sep = '\t',
                            data.table = FALSE, fill = TRUE, check.names=FALSE,
                            na.string=c("NA"), nThread = 10)
      head(df)
      toKeep<-colnames(df)[2]
      rnaCounts<-cbind(rnaCounts,df[,toKeep])
      toKeep<-unlist(strsplit(toKeep, '_', fixed=TRUE))[1]
      colNames<-append(colNames,toKeep)
  }
  
  head(rnaCounts)
  colNames<-gsub("TS-677","S-lycLYC3153",colNames)
  colNames<-gsub("TS-684","S-lycEA02054",colNames)
  colNames<-gsub("TS-731","S-lycEA00940",colNames)
  tmp<-grep("^P", colNames, value = T)
  tmp<-gsub('-1','', tmp)
  tmp2<-grep("^P", colNames, value = T, invert = T)
  colNames<-c(tmp, tmp2)
  batch<-data.table::fread("batch_separation.csv", sep = '\t',
                             data.table = FALSE, fill = TRUE, check.names=FALSE,
                             na.string=c("NA"), nThread = 10)
  rownames(batch)<-batch$V1
  setdiff(colNames, batch$V1)
  batch<-batch[colNames,]
    
  ##### adjust by batch effects
  adjustedCounts <- ComBat_seq(rnaCounts, batch=batch$V2, group=NULL)
  adjustedCounts<-cbind.data.frame(adjustedCounts,df[,1])
  ####### re name some samples
  colnames(adjustedCounts)<-c(colNames,'ID')
    
  refSample<-grep('TS-253', colNames, value = T)
  allSample<-grep('TS-253', colNames, value = T, invert = T)
  rownames(adjustedCounts)<-adjustedCounts$ID
  batch<-batch[c(refSample,allSample),]
    
  adjustedCounts<-adjustedCounts[,c(refSample,allSample)]
  myCondition<-data.frame (SampleID  = c(refSample, allSample),
                             Condition = c('Control', rep('Treat', times=length(allSample))))
    # # do the deseq transformation 
    #adjusted
  dds_adjusted <- DESeqDataSetFromMatrix(countData = adjustedCounts, 
                                           colData = myCondition, design = ~Condition)
  dds_adjusted <- DESeq(dds_adjusted)
  res_adjusted<-results(dds_adjusted)
  ########################
  ## Normalization done with DESeq
  df<-as.data.frame(counts(dds_adjusted, normalized=TRUE))
  df$Geneid<-rownames(dds_adjusted)
  samples2keep<-grep('Geneid', colnames(df), value = T, invert = T)
  df$Geneid <- sub("\\.\\d+$", "", df$Geneid)
  head(df)
  data.table::fwrite(df, file= 'leaf_transcriptome-counts_normalized_batch-adjusted.tsv', quote=F,
                       row.names=F,col.names = T,sep="\t")
