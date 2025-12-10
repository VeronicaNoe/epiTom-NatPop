suppressPackageStartupMessages({
  require("data.table")
  library(ggplot2)
  library(dplyr)
})
setwd("/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation")
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"

############
input<-list.files(pattern = "^ch", full.names = FALSE)
out<-c()
for (i in input){
  dmr<-unlist(strsplit(i,".", fixed=T))[1]
  dmr<-unlist(strsplit(dmr,"_", fixed=T))[2]
  df<- data.table::fread(i, 
                        sep = '\t', data.table = T, header = T,
                        fill = TRUE, na.string=c("NA"), nThread = 20) 
  df$V204<-NULL
  sample2keepCol<-grep('china',colnames(df), invert = T, value = T )
  
  df<-df %>%
    select(infoCol) %>%
    dplyr::mutate(nSamples = rowSums(!is.na(df[,..sample2keepCol])))
  df$dmr<-dmr
  out<-rbind.data.frame(out, df)
}

window_size <- c(0.5,1.0,1.5,2,2.5)  # Define your window size
out$start<-out$start/1000000
out$chr<-gsub('SL2.50ch','', out$chr)

densityByBins <- out %>%
  group_by(dmr) %>%
  mutate(totalDMRs = n()) %>%
  ungroup() %>%
  mutate(bin = floor(start / 0.1) * 0.1) %>%
  group_by(chr, bin, dmr, totalDMRs) %>%
  summarise(nDMR = n(), .groups = 'drop') %>%
  mutate(dmrProportion = nDMR / totalDMRs)

ggplot(densityByBins, aes(x = bin, y = dmrProportion, fill = dmr, color=dmr)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~chr, scales = "free_x") +
  scale_color_manual(values=c("#820a86", "#ffca7b"))+
  scale_fill_manual(values=c("#820a86", "#ffca7b"))+
  labs(x = "Bin", y = "Density of DMRs", fill = "DMR Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()  )
ggsave(paste0(outPlot,"03.06_DMR_density_allChr.pdf"), width = 60, height = 40, units = "cm")

tmp<-densityByBins%>%
  filter(chr=="01")
ggplot(tmp, aes(x = bin, y = dmrProportion, fill = dmr, color=dmr)) +
  geom_bar(stat = "identity", position = "dodge") +
  #facet_wrap(~chr, scales = "free_x") +
  scale_color_manual(values=c("#820a86", "#ffca7b"))+
  scale_fill_manual(values=c("#820a86", "#ffca7b"))+
  labs(x = "Bin", y = "Density of DMRs", fill = "DMR Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()  )
ggsave(paste0(outPlot,"03.07_DMR_density_chr01.pdf"), width = 60, height = 40, units = "cm")
