suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggridges)
  
})
#setwd("/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation")
setwd("/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection")
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"

input<-list.files(pattern = "DMR.methylation", full.names = FALSE)
toKeep<-c("TE-wo-gene", 'gene')
input<- input[grep(paste(toKeep, collapse = "|"), input)]

for(k in toKeep){
  if(k =="gene"){
    inData<-grep(k, input, value = T)
    inData<-grep("TE-wo", inData, value = T, invert = T)
  } else {
    inData<-grep(k, input, value = T)
    }
  out<-c()
  for(i in inData){
    toAnno<-gsub(".methylation", "",i)
    DMR<-unlist(strsplit(toAnno, ".", fixed = TRUE))[2]
    anno<-unlist(strsplit(toAnno, ".", fixed = TRUE))[1]
    df<- data.table::fread(i, sep = '\t', data.table = TRUE, fill = TRUE, na.string=c("NA"), nThread = 20)
    df <- df %>%
      pivot_longer(
        cols = -c(V1,V2,V3),
        names_to = "samples",
        values_to = "value")
    gc()
    colnames(df)<-c("chr",'start','end','sample','methLevel')
    if(k =="gene"){
      anno<-k
    }
    df$anno<-anno  
    df$DMR<-DMR
    out<-rbind.data.frame(out, df)
    gc()
  }
  p<-ggplot(out, aes(x = methLevel, color=DMR, fill=DMR)) + 
    geom_density(alpha=0.4) +
    labs(title="",x="Methylation level", y = "Density")+
    scale_fill_manual(values=c("C-DMR"="#820a86", "CG-DMR"="#ffca7b"))+
    scale_color_manual(values=c("C-DMR"="#820A86", "CG-DMR"="#ffca7b"))+
    theme_ridges(font_size = 8, center_axis_labels = T)+
    theme_minimal()
  ggsave(paste0(outPlot,"05.05_",anno,"_C_vs_CG_DMR_methLevels_density.pdf"),
         plot=p,width = 60, height = 60, units = "cm")
  gc()
}
