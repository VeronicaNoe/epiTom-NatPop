suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
#setwd("/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/data")
#outPlot<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/plots/"
#outDir<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/results/"
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig1/ab_global-methylation")
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"

############################################################
##### methylation per sample
#################################################################
input<-list.files(pattern = ".windowed", full.names = FALSE)
column_names <- c('acc','tissue','context','methylation')
out <- data.frame(matrix(ncol = length(column_names), nrow = length(input)))
colnames(out)<-column_names
for (i in 1:length(input)){
  acc<-unlist(strsplit(input[i],"_", fixed=T))[1]
  tissue<-unlist(strsplit(input[i],"_", fixed=T))[2]
  context<-gsub('.windowed','', unlist(strsplit(input[i],"_", fixed=T))[4])
  df<- data.table::fread(paste0(acc,"_",tissue,"_Biseq_",context,".windowed"), 
                           sep = '\t', data.table = FALSE, header = F,
                           fill = TRUE, na.string=c("NA"), nThread = 20) 
  
  colnames(df)<-c('chr','start','end','meth', 'unmeth')
  totalMeth<-sum(df$meth)
  totalUnMeth<-sum(df$unmeth)
  total<-sum(totalMeth,totalUnMeth)
  out[i,'acc']<-acc
  out[i,'tissue']<-tissue
  out[i,'context']<-context
  out[i,'methylation']<-totalMeth/total
}
out_wide <- out %>%
  pivot_wider(names_from = context, values_from = methylation)
write.table(out_wide, paste0(outDir, "02.01_global_methylation_sample.tsv"), sep = '\t', row.names = FALSE, quote = FALSE)

toPlot<-subset(out, tissue=='leaf')
toPlot$acc<-gsub('S-penLA0716', 'S-penLA1272',toPlot$acc)
toPlot$acc<-gsub('LA2838A', 'S-lycLA2838A',toPlot$acc)

toRemove<-c('LA0534','Ohio-8245','S-lycLA2845','S-pimLA1578','TS-244')
toPlot<-toPlot[!(toPlot$acc %in% toRemove),]

mean_values <- toPlot %>%
  group_by(context, tissue) %>%
  dplyr::summarise(mean_methylation = mean(methylation))
mean_values
# # A tibble: 3 × 3
# # Groups:   context [3]
# context tissue mean_methylation
# <chr>   <chr>             <dbl>
# 1 CG      leaf             0.643 
# 2 CHG     leaf             0.414 
# 3 CHH     leaf             0.0642
summary_stats <- toPlot %>%
  group_by(context) %>%
  summarise(
    mean_methylation = mean(methylation, na.rm = TRUE),
    median_methylation = median(methylation, na.rm = TRUE),
    sd_methylation = sd(methylation, na.rm = TRUE),
    min_methylation = min(methylation, na.rm = TRUE),
    max_methylation = max(methylation, na.rm = TRUE),
    cv_methylation = (sd(methylation, na.rm = TRUE) / mean(methylation, na.rm = TRUE)) * 100
  )
print(summary_stats, width=Inf)
# # A tibble: 3 × 7
# context mean_methylation median_methylation sd_methylation min_methylation
# <chr>              <dbl>              <dbl>          <dbl>           <dbl>
# 1 CG                0.643              0.649          0.0582          0.439 
# 2 CHG               0.414              0.409          0.0515          0.252 
# 3 CHH               0.0642             0.0640         0.0137          0.0376
# max_methylation cv_methylation
# <dbl>          <dbl>
# 1           0.796         9.05
# 2           0.610         12.4 
# 3           0.104         21.4 

treeOrder<-data.table::fread(paste0(outDir,"01.01_data_general_SNPs_order.csv"), 
                             sep = '\t', data.table = FALSE, header = T,
                             fill = TRUE, na.string=c("NA"), nThread = 20) 
treeOrder$label<-gsub('S_','S-',treeOrder$label)
treeOrder$label<-gsub('_1','',treeOrder$label)
treeOrder$label<-gsub('BGV006775','S-pimBGV006775',treeOrder$label)

toPlot$orderAcc <- factor(toPlot$acc, levels = treeOrder$label)

ggplot(toPlot, aes(x = orderAcc, y = methylation, color = context, fill=context)) +
  geom_bar(stat="identity", width=0.5, alpha = 0.4) +
  labs(x = "Analyzed accessions", 
       y = "Methylation levels")+
  scale_color_manual(values=c("#CCCC00","#0033CC", "#CC0066"))+
  scale_fill_manual(values=c("#CCCC00","#0033CC", "#CC0066"))+
  facet_grid(context~tissue, scales='free')+
  #ylim(0,1)+
  geom_hline(data = mean_values, aes(yintercept = mean_methylation), color = "black",
             linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),panel.background = element_blank())
ggsave(paste0(outPlot,"02.01.Accession_mean-global-methylation.pdf"),  width = 60, height = 20, units = "cm")  

# same with points
# ggplot(toPlot, aes(x = orderAcc, y = methylation, color = context, fill=context)) +
#   geom_point(size = 4, alpha = 0.6,shape = 21, stroke = 0.5) +  # Use geom_point for points
#   labs(x = "Analyzed accessions",
#        y = "Methylation levels")+
#   scale_color_manual(values=c("black","black", "black"))+
#   scale_fill_manual(values=c("#CCCC00","#0033CC", "#CC0066"))+
#   facet_grid(context~tissue, scales='free')+
#   #ylim(0,1)+
#   geom_hline(data = mean_values, aes(yintercept = mean_methylation), color = "black",
#              linetype = "dashed")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),panel.background = element_blank())
# ggsave(paste0(outPlot,"02.02.Accession_mean-global-methylation_points.pdf"),  width = 60, height = 20, units = "cm")
