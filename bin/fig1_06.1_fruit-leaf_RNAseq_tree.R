suppressPackageStartupMessages({
  library(data.table)
  library(sva) #combat
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(vegan)
  library(dendextend)
})
setwd("/mnt/disk2/vibanez/07_rnaseq-processing/07.3_counts")
outDir<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/results/"
outPlot<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/plots/"

############################################################
############### load expression count
############################################################
tissue<-c("leaf","fruit")
rnaCounts <- c()
geneNames<-c()
colNames<-c()
for( t in tissue){
  inDir<-paste0(t,"/")
  input<-list.files(path=inDir, pattern = ".tsv", full.names = FALSE)
  for (i in input){
    sampleName<-gsub('_RNAseq.counts_edited.tsv','',i)
    df <- data.table::fread(paste0(inDir, i), sep = '\t', skip = "Geneid",
                            data.table = T,fill = TRUE, check.names = FALSE, 
                            na.strings = c("NA"), nThread = 10)
    df<-df %>%
      filter(!grepl("Solyc00g", Geneid) & !grepl("TomatoPan", Geneid))
    colnames(df)[2]<-sampleName
    rnaCounts<-cbind(rnaCounts,df[,..sampleName])
    geneNames<-df$Geneid
    colNames<-append(colNames,sampleName)
  }
}
head(rnaCounts)
colNames<-gsub("TS-677_leaf","S-lycLYC3153_leaf",colNames)
colNames<-gsub("TS-684_leaf","S-lycEA02054_leaf",colNames)
colNames<-gsub("TS-731_leaf","S-lycEA00940_leaf",colNames)

colnames(rnaCounts)<-colNames
# Adjust leaf counts by batch 
leaf_samples<-grep('leaf', colNames, value = T)

batch<-data.table::fread("batch_separation.csv", sep = '\t',
                         data.table = FALSE, fill = TRUE, check.names=FALSE,
                         na.string=c("NA"), nThread = 10)
rownames(batch)<-batch$V1
setdiff(leaf_samples, batch$V1)
batch<-batch[leaf_samples,]

##### adjust by batch effects
leaf_data<-rnaCounts[,..leaf_samples]
adjustedCounts <- ComBat_seq(leaf_data, batch=batch$V2, group=NULL)
adjustedCounts<-cbind.data.frame(adjustedCounts,df[,1])

####### 
fruit_samples<-grep('fruit', colNames, value = T)
counts<-cbind(rnaCounts[,..fruit_samples],adjustedCounts)

refSample<-'TS-253_leaf'
allSample<-grep(refSample, conames(counts), value = T, invert = T)

counts<-counts[,c(refSample,allSample)]
rownames(counts)<-counts$ID

myCondition<-data.frame (SampleID  = c(refSample, allSample),
                         Condition = c('Control', rep('Treat', times=length(allSample))))
# # do the deseq transformation 
#adjusted
dds_adjusted <- DESeqDataSetFromMatrix(countData = counts, 
                                       colData = myCondition, design = ~Condition)
dds_adjusted <- DESeq(dds_adjusted)
res_adjusted<-results(dds_adjusted)
########################
## Normalization done with DESeq
#unadjusted
df<-as.data.frame(counts(dds_adjusted, normalized=TRUE))
df$Geneid<-rownames(dds_adjusted)
samples2keep<-grep('Geneid', colnames(df), value = T, invert = T)
df$Geneid <- sub("\\.\\d+$", "", df$Geneid)
rownames(df)<-df$Geneid
head(df)
###############################

  
  samples2keep<-c("S-lycEA00375_leaf", "S-lycLYC2962_leaf","S-lycPI129097_leaf","S-lycPI311117_leaf",
                  "S-lycSG16_leaf", "TS-154_leaf", "TS-20_leaf", "TS-21_leaf",
                  "TS-223_leaf", "TS-238_leaf", "TS-253_leaf", "TS-528_leaf", "TS-533_leaf", 
                  "TS-6_leaf", "S-lycEA00375_fruit", "S-lycLYC2962_fruit",
                  "S-lycPI129097_fruit","S-lycPI311117_fruit", "S-lycSG16_fruit", 
                  "TS-154_fruit", "TS-20_fruit", "TS-21_fruit","TS-223_fruit", "TS-238_fruit", 
                  "TS-253_fruit", "TS-528_fruit", "TS-533_fruit", "TS-6_fruit")
  # 1. Subset rnaCounts using samples2keep, and reorder based on refSample and allSample
  rnaCounts <- rnaCounts[, ..samples2keep]
  refSample <- grep('TS-253_leaf', samples2keep, value = TRUE)
  allSample <- grep('TS-253_leaf', samples2keep, value = TRUE, invert = TRUE)
  
  # Reorder columns: refSample first, then allSample
  new_order <- c(refSample, allSample)
  rnaCounts <- rnaCounts[, new_order, with = FALSE]
  
  # Convert rnaCounts to data frame and set gene names as row names
  setDF(rnaCounts)
  rownames(rnaCounts) <- geneNames
  # 2. Define conditions for each sample
  myCondition <- data.frame(SampleID  = c(refSample, allSample),
    Condition = c('Control', rep('Treat', times = length(allSample))) )
  
  # 3. Create DESeq2 dataset
  dds_adjusted <- DESeqDataSetFromMatrix(
    countData = rnaCounts,
    colData = myCondition,
    design = ~ Condition
  )
  # 4. Run DESeq normalization and differential analysis
  dds_adjusted <- DESeq(dds_adjusted)
  # 5. Extract normalized counts
  normalized_counts <- counts(dds_adjusted, normalized = TRUE)
  # 6. (Optional) If you want the results of differential expression analysis
  res_adjusted <- results(dds_adjusted)
  # 7. View normalized counts
  head(normalized_counts)
  # (Optional) Convert normalized counts to a data frame
  df_normalized <- as.data.frame(normalized_counts)
  
  ### do the dendro
  sampleTab<-data.table::fread("/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/data/00_sampleTab.tsv", header = TRUE,sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
  head(sampleTab)
  toGroup<-sampleTab[sampleTab$MethName %in% samples2keep,]
  rownames(toGroup)<-toGroup$MethName
  # Define your group and color mapping
  toGroup <- toGroup %>%
    mutate(Col = case_when(Group == "WILD" ~ "#000000ff",
                           Group == "PIM"  ~ "#00a681ff",
                           Group == "SLC"  ~ "#ed83b5ff",
                           Group == "SLL"  ~ "#00a4deff",
                           TRUE            ~ "grey" ),
           tissueCol = case_when(Organ == "leaf" ~ "#28AC34",
                           Organ == "fruit"  ~ "#cf4f62"))  
  
    toDendro<-t(df_normalized)
    rownames(toDendro)<-samples2keep
    toDendro[1:20,1:10]
    #get distance
    dis <-vegdist(toDendro, na.rm=TRUE, "euclid")
    #get dendro
    dend<-as.dendrogram(hclust(as.dist(dis), "ward.D2"))
    # get same order
    toOrder<-dend %>% labels # get the labels of the tree to sort
    toGroup<-toGroup[toOrder,]
    toGroup <- toGroup[match(toOrder, toGroup$MethName), ]
    pdf(paste0(outPlot,"leaf-fruit_dendrogram_RNAseq.pdf"),width = 15, height = 8)
    dend %>% 
      set("labels_col", toGroup$Col) %>%  
      set("labels_cex", 0.5) %>%
      set("leaves_pch", 19) %>% 
      set("leaves_cex", 0.5) %>% 
      set("leaves_col", toGroup$tissueCol) %>%
      plot()
    dev.off()
    