# load tidyverse package
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# setwd("/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/data")
# outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
# outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig1/ac_pca")
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"


input<-list.files(pattern = "eigenval", full.names = FALSE)
### PCA
for(i in 1:length(input)){
  sample<-gsub('.eigenval','',input[i])
  anno<-unlist(strsplit(sample, "_", fixed = TRUE))[1]
  marker<-unlist(strsplit(sample, "_", fixed = TRUE))[2]
  sampleTab<-data.table::fread("/mnt/disk2/vibanez/10_data-analysis/Fig1/00_sampleTab.tsv", sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 16) 
  sampleTab<-subset(sampleTab, Organ=="leaf")
  if(marker == 'SNPs'){
    sampleTab <- subset(sampleTab, SNPnames != "")
    sampleLabel<-unique(sampleTab$SNPnames)
    sampleTab<-sampleTab[sampleTab$SNPnames %in% sampleLabel,]
    rownames(sampleTab)<-sampleLabel
  }else{
    rownames(sampleTab)<-sampleTab$MethName
    group<-sampleTab[sampleTab$MethName,"Group"]
  }
  # read in data
  pca <- read.table(paste0(sample,".eigenvec"))
  eigenval <- scan(paste0(sample,".eigenval"))
  # sort out the pca data
  # remove nuisance column
  if(marker == 'SNPs'){
    pca$V1 <-NULL 
    # set names
    names(pca)[1] <- "ind"
    names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
    group<-sampleTab[pca$ind,"Group"]
    pca$Group<-group
  }else{
    pca$V2 <-NULL 
    names(pca)[1] <- "ind"
    names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
    group<-sampleTab[paste0(pca$ind,'_leaf'),"Group"]
    pca$Group<-group
  }
  
  # sort out the individual species and pops
  # first convert to percentage variance explained
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  #define groups
  
  # # calculate the cumulative sum of the percentage variance explained
  cumsum(pve$pve)
  rownames(pca)<-pca$ind
  pca$Group <- factor(pca$Group, levels = c("WILD", "PIM", "SLC", "SLL"))
  # plot pca
  ggplot(pca, aes(PC1, PC2, label = ind, color=Group)) +
    geom_point(size = 5)+
    scale_color_manual(values=c("#000000ff","#00a681ff","#ed83b5ff","#00a4deff"))+
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
    ggtitle(paste0(marker," ",anno))+
    #geom_text(aes(label=ind), hjust=0, vjust=1)+  
    theme(plot.title = element_text(size = 10))  # Adjust font size here
  ggsave(paste0(outPlot,"04.01.",marker,"_",anno,"_PCA.pdf"),  width = 20, height = 20, units = "cm")
}
