########## merge data
#do this in the Camus
suppressPackageStartupMessages({
  require("data.table")
  library(dendextend)
  library(dplyr)
  library(vegan)
  library(tidyr)
  library(ggplot2)
})
  
setwd("/mnt/disk2/vibanez/08_fruit-processing/ai_tissue-merge/result")
input<-list.files(pattern = ".bed", full.names = FALSE)

dendro<-c()
for(i in 1:length(input)){
  df<- data.table::fread(input[i], sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
  sample<-unlist(strsplit(input[i],".", fixed = TRUE))[1]
  dmr<-unlist(strsplit(sample,"_", fixed = TRUE))[2]
  colnames(df)<-gsub("_Biseq","",colnames(df))
  #create each vector
  tmp<-cbind.data.frame(dmr,df)
  dendro<-rbind.data.frame(dendro, tmp)
}
head(dendro)
colnames(dendro) <- gsub("S-lycEA00375_leaf$", "TS-528_leaf_R", colnames(dendro))
colnames(dendro) <- gsub("^TS-528_leaf$", "S-lycEA00375_leaf_R", colnames(dendro))
colnames(dendro) <- gsub("_R$", "", colnames(dendro))

samples2keep<-c("S-lycEA00375_leaf", "S-lycLYC2962_leaf","S-lycPI129097_leaf","S-lycPI311117_leaf",
                "S-lycSG16_leaf", "TS-154_leaf", "TS-20_leaf", "TS-21_leaf",
                "TS-223_leaf", "TS-238_leaf", "TS-253_leaf", "TS-528_leaf", "TS-533_leaf", 
                "TS-6_leaf", "S-lycEA00375_fruit", "S-lycLYC2962_fruit",
                "S-lycPI129097_fruit","S-lycPI311117_fruit", "S-lycSG16_fruit", 
                "TS-154_fruit", "TS-20_fruit", "TS-21_fruit","TS-223_fruit", "TS-238_fruit", 
                "TS-253_fruit", "TS-528_fruit", "TS-533_fruit", "TS-6_fruit")
dendro$V215<-NULL
# data.table::fwrite(dendro, file= "DMRs_leaf-fruit.bed", quote=F,
#                    row.names=F,col.names = T,sep="\t")
infoCol<-c("dmr", "chr", "start", "end")
toKeep<-c(infoCol, samples2keep)
sampleTab<-data.table::fread("00_sampleTab.tsv", header = TRUE,sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
head(sampleTab)
toGroup<-sampleTab[sampleTab$Sample2 %in% samples2keep,]
rownames(toGroup)<-toGroup$Sample2
# Define your group and color mapping
toGroup <- toGroup %>%
  mutate(Col = case_when(
    Group == "WILD" ~ "#000000ff",
    Group == "PIM"  ~ "#00a681ff",
    Group == "SLC"  ~ "#ed83b5ff",
    Group == "SLL"  ~ "#00a4deff",
    TRUE            ~ Col  # Keep the original color if no match
  ))  

## subset for DMR
DMR<-c("C-DMR", "CG-DMR")
for(d in DMR){
  toDendro<-subset(dendro, dmr==d)
  toMeth<-toDendro[,toKeep]
  
  average_methylation <- toMeth %>%
    group_by(dmr) %>%
    summarise(across(all_of(samples2keep), mean, na.rm = TRUE))
  
  average_methylation <- average_methylation %>%
    pivot_longer(cols = -dmr,   # Keep 'dmr' as the fixed column
                 names_to = "Sample",    # The column where sample names go
                 values_to = "MethyLevels") %>%
    mutate(Sample_ID = sub("_.*", "", Sample),   # Extract the Sample_ID
           Tissue = sub(".*_", "", Sample))
  
  toDendro<-toDendro[,samples2keep]
  toDendro<-t(toDendro)
  rownames(toDendro)<-samples2keep
  #toDendro[1:20,1:10]
  #get distance
  dis <-vegdist(toDendro, na.rm=TRUE, "euclid")
  #get dendro
  dend<-as.dendrogram(hclust(as.dist(dis), "ward.D2"))
  # get same order
  toOrder<-dend %>% labels # get the labels of the tree to sort
  toGroup<-toGroup[toOrder,]
  toGroup <- toGroup[match(toOrder, toGroup$Sample2), ]
  pdf(paste0(d,"_leaf-fruit_dendrogram.pdf"),width = 15, height = 8)
  dend %>% set("labels_col", toGroup$Col) %>%  set("labels_cex", 0.5) %>%
    set("leaves_pch", 19) %>% set("leaves_cex", 0.5) %>% set("leaves_col", toGroup$ColTissue) %>%
    plot()
  dev.off()
  
  average_methylation<-as.data.frame(average_methylation)
  rownames(average_methylation)<-average_methylation$Sample
  average_methylation<-average_methylation[toOrder,]
  average_methylation$Sample <- factor(average_methylation$Sample, levels = rownames(average_methylation))
  
  ggplot(average_methylation, aes(x = Sample, y = MethyLevels, color = Tissue, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5, alpha = 0.4) +
    labs(x = "Analyzed accessions", 
         y = "Methylation levels") +
    scale_color_manual(values = c("leaf" = "#28AC34", "fruit" = "#cf4f62")) +
    scale_fill_manual(values = c("leaf" = "#28AC34", "fruit" = "#cf4f62")) +
    facet_grid(dmr ~ ., scales = 'free') +
    ylim(0, 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank())
  ggsave(paste0(d,"_accession-tissue_mean-global-methylation.pdf"), width = 60, height = 20, units = "cm")
}


