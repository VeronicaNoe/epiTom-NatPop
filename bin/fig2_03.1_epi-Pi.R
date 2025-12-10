#!/users/bioinfo/vibanez/miniconda3/envs/methylKit/bin/Rscript
suppressPackageStartupMessages({
  require("data.table")
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(ggpubr)
  library(emmeans)
})
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig2/ac_calculate-nucleotide-dmr-diversity/bc_merge-column-wise-group/")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/plots/"

infoCols<-c("chr", "start", "end")
groupNames<-c('1-WILD','2-PIM','3-SLC','4-SLL')
#############################################
## A get overall diversity index
#############################################
# boxplots
{
anno<-c("gene","gene-TE","intergenic","promotor", "TE-wo-gene")
diversity<-c()
for(a in anno){
  pi<- data.table::fread(paste0('SNPs_',a,".merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  #head(pi)
  pi$V8<-NULL
  colnames(pi)<-c(infoCols,'1-WILD_SNPs','2-PIM_SNPs','3-SLC_SNPs','4-SLL_SNPs')
  cDMR<- data.table::fread(paste0('C-DMR_',a,"_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  cDMR$V8<-NULL
  colnames(cDMR)<-c(infoCols,'1-WILD_C-DMR','2-PIM_C-DMR','3-SLC_C-DMR','4-SLL_C-DMR')
  cgDMR<- data.table::fread(paste0('CG-DMR_',a,"_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  cgDMR$V8<-NULL
  colnames(cgDMR)<-c(infoCols,'1-WILD_CG-DMR','2-PIM_CG-DMR','3-SLC_CG-DMR','4-SLL_CG-DMR')
  
  merged_table <- dplyr::left_join(pi, cDMR, by = infoCols)
  merged_table <- dplyr::left_join(merged_table, cgDMR, by = infoCols)
  head(merged_table)
  merged_table$annotation<-a
  diversity<-rbind.data.frame(diversity, merged_table)
}
head(diversity)
diversity <- reshape2::melt(diversity, id.vars = c("chr", "start", "end", "annotation"), variable.name = "samples",
                    value.name = "diversity")
diversity <- tidyr::separate(diversity, samples, into = c("group", "marker"), sep = "_", remove = T)
diversity$diversity<-diversity$diversity*1000
data.table::fwrite(diversity, file=paste0(outDir,"03.1_diversity_merged.tsv"), 
                   quote=F,row.names=F,col.names = T,sep="\t")
}
# diversity<- data.table::fread(paste0(outDir,"diversity_merged.tsv"), 
#                               sep = '\t', data.table = FALSE, fill = TRUE, 
#                               na.string=c("NA"), nThread = 20)
quantiles <- quantile(diversity$diversity, probs = seq(0, 1, 0.01), na.rm = TRUE)
#diversity$diversity[diversity$diversity > quantiles[96]] <- NA # missing for Q90= quantiles[91] | 90% |22.5907 
diversity <- na.omit(diversity)
head(diversity)

# boxplot for diversity
{
ggplot(diversity, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
  #geom_violin( na.rm = T) +
  geom_boxplot(width=0.5, na.rm = T, 
               alpha=0.1,lwd = 0.3,outlier.shape = NA)+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) + 
  labs(x = "", 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_fill_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  scale_color_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~annotation, scales = 'free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))+
  theme_bw()
ggsave(paste0(outPlot,"03.01.diversity-per-features.pdf"),  width = 21, height = 20, units = "cm")  

ggplot(diversity, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
  #geom_violin( na.rm = T) +
  geom_boxplot(width=0.5, na.rm = T, 
               alpha=0.1,lwd = 0.3,outlier.shape = NA)+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) + 
  labs(x = "", 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_fill_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  scale_color_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~., scales = 'free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))+
  theme_bw()
ggsave(paste0(outPlot,"03.02.diversity-all-features.pdf"),  width = 21, height = 20, units = "cm")  

tmp<-subset(diversity, group!='1-WILD')
# summary_stats <- tmp %>%
#   group_by(group, annotation, marker) %>%
#   summarise(mean_diversity = mean(diversity),
#             sd_diversity = sd(diversity),
#             max_diversity = max(diversity),
#             min_diversity = min(diversity),
#             n_diversity = n())
# testAnova <- aov(diversity ~ group + annotation + marker, data = tmp)
# testTukey <- TukeyHSD(testAnova)
# glm(diversity ~ group, data = tmp1, family = gaussian)
# # View the summary statistics
# summary_stats
#cmpr <- list(c("3-SLC","2-PIM"),c("3-SLC","4-SLL"),c("4-SLL","2-PIM"))
ggplot(tmp, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
  #geom_violin(na.rm = T) +
  geom_boxplot(width=0.5,alpha=0.1,na.rm = T,lwd = 0.3,outlier.shape = NA)+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) + 
  #stat_compare_means(comparisons = cmpr)+ # Add pairwise comparisons p-value
  #stat_compare_means()  +   # Add global p-value
  labs(x = "", 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_fill_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  scale_color_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~annotation, scales = 'free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))+
  theme_bw()
ggsave(paste0(outPlot,"03.03.diversity-per-features_woWILD.pdf"),  width = 21, height = 20, units = "cm")  

ggplot(tmp, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
  #geom_violin(na.rm = T) +
  geom_boxplot(width=0.5,alpha=0.1,na.rm = T,lwd = 0.3,outlier.shape = NA)+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) + 
  #stat_compare_means(comparisons = cmpr)+ # Add pairwise comparisons p-value
  #stat_compare_means()  +   # Add global p-value
  labs(x = "", 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_fill_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  scale_color_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~., scales = 'free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))+
  theme_bw()
ggsave(paste0(outPlot,"03.04.diversity-all-features_woWILD.pdf"),  width = 21, height = 20, units = "cm")  

}
#############################################
## B plot diversity along the chromosomes
#############################################
pi2plot<-subset(diversity, chr!='SL2.50ch00')
pi2plot$chr<-gsub('SL2.50','',pi2plot$chr)
pi2plot$start<-pi2plot$start/1000000
pi2plot$end<-NULL
pi2plot$annotation<-NULL
head(pi2plot)
# plot along chr
window_size <- 2.5  # Define your window size
pi2plot <- pi2plot %>%
  group_by(chr, group, marker, window = floor(start / window_size)) %>%
  summarize(mean_diversity = mean(diversity, na.rm = TRUE),
            start = first(start)) %>%
  ungroup()
# Plot the mean diversity by windows
ggplot(pi2plot, aes(x = start, y = mean_diversity, group = group)) +
  geom_line(aes(color = group), size = 0.5) +
  labs(x = "chromosome position (Mb)", 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_color_manual(values=c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~chr, scales='free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw()
ggsave(paste0(outPlot,"03.05.diversity_along-chr.pdf"),  width = 61, height = 20, units = "cm")  
#woWILD
tmp<-subset(pi2plot, group!='1-WILD')
ggplot(tmp, aes(x = start, y = mean_diversity, group = group)) +
  geom_line(aes(color = group), size = 0.5) +
  labs(x = "chromosome position (Mb)", 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_color_manual(values=c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~chr, scales='free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw()
ggsave(paste0(outPlot,"03.06.diversity_along-chr_woWILD.pdf"),  width = 61, height = 20, units = "cm")  

# ANOVA
tmp1<-subset(tmp, marker=='C-DMR')
glm_model <- glm(diversity ~ group, data = tmp1, family = gaussian)
pairwise_result <- emmeans(glm_model, pairwise ~ group)
### NOT NECESSARU
{# pi<- data.table::fread(paste0("SNPs_general_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
# #head(pi)
# pi$V8<-NULL
# colnames(pi)<-c(infoCols,'1-WILD_SNPs','2-PIM_SNPs','3-SLC_SNPs','4-SLL_SNPs')
# cDMR<- data.table::fread(paste0("C-DMR_general_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
# cDMR$V8<-NULL
# colnames(cDMR)<-c(infoCols,'1-WILD_C-DMR','2-PIM_C-DMR','3-SLC_C-DMR','4-SLL_C-DMR')
# cgDMR<- data.table::fread(paste0("CG-DMR_general_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
# cgDMR$V8<-NULL
# colnames(cgDMR)<-c(infoCols,'1-WILD_CG-DMR','2-PIM_CG-DMR','3-SLC_CG-DMR','4-SLL_CG-DMR')
# 
# allFeaturesDiversity <- dplyr::left_join(pi, cDMR, by = infoCols)
# allFeaturesDiversity <- dplyr::left_join(allFeaturesDiversity, cgDMR, by = infoCols)
# 
# allFeaturesDiversity <- reshape2::melt(allFeaturesDiversity, id.vars = c("chr", "start", "end"), variable.name = "samples",
#                             value.name = "diversity")
# allFeaturesDiversity <- tidyr::separate(allFeaturesDiversity, samples, into = c("group", "marker"), sep = "_", remove = T)
# allFeaturesDiversity$diversity<-allFeaturesDiversity$diversity*1000
# 
# quantiles <- quantile(allFeaturesDiversity$diversity, probs = seq(0, 1, 0.01), na.rm = TRUE)
# #allFeaturesDiversity$diversity[allFeaturesDiversity$diversity > quantiles[96]] <- NA # missing for Q90= quantiles[91] | 90% |22.5907 
# #allFeaturesDiversity <- na.omit(allFeaturesDiversity)
# 
# ggplot(allFeaturesDiversity, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
#   #geom_violin( na.rm = T) +
#   geom_boxplot(width=0.5, na.rm = T,
#                alpha=0.1,lwd = 0.3,outlier.shape = NA)+
#   stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) +
#   labs(x = "", 
#        y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
#   scale_fill_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
#   scale_color_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
#   facet_grid(marker~., scales = 'free_y')+
#   theme_bw() +
#   theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))
# ggsave(paste0(outPlot,"diversity_general.pdf"),  width = 20, height = 20, units = "cm")
# 
# tmp<-subset(allFeaturesDiversity, group!="1-WILD")
# ggplot(tmp, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
#   #geom_violin( ) +
#   geom_boxplot(width=0.5, 
#                alpha=0.1,lwd = 0.3,outlier.shape = NA)+
#   stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) +
#   #stat_compare_means(comparisons = cmpr)+ # Add pairwise comparisons p-value
#   #stat_compare_means()  +   # Add global p-value
#   labs(x = "", 
#        y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
#   scale_fill_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
#   scale_color_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
#   facet_grid(marker~., scales = 'free_y')+
#   theme_bw() +
#   theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))
# ggsave(paste0(outPlot,"diversity_general_woWild.pdf"),  width = 20, height = 20, units = "cm")
}

{
#############################################
## C plot diversity fro selective sweeps
#############################################
DS<-subset(diversity, chr == "SL2.50ch12")
DS<-subset(DS, group != "1-WILD")
DS$start<-DS$start/1000000
DS<-subset(DS, start >= 0 )
DS<-subset(DS, start <= 10)
DS <- DS %>%
  group_by(chr, group,marker,window = floor(start / 0.1)) %>%
  summarize(mean_diversity = mean(diversity, na.rm = TRUE),
            start = first(start)) %>%
  ungroup()
ggplot(DS, aes(x = start, y = mean_diversity, group = group)) +
  geom_line(aes(color = group), size = 0.5) +
  labs(x = paste0("chromosome 12 position (Mb)"), 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_color_manual(values=c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~., scales='free')+
  scale_x_continuous(breaks =  seq(0, 10, by = 1))+
  theme_bw()
ggsave(paste0(outPlot,"03.07.diversity_along-DomesticatedSweep.pdf"),  width = 21, height = 20, units = "cm")  

IS<-subset(diversity, chr == "SL2.50ch02")
IS<-subset(IS, group != "1-WILD")
IS$start<-IS$start/1000000
IS<-subset(IS, start >= 34 )
IS<-subset(IS, start <= 50)
IS <- IS %>%
  group_by(chr, group,marker,window = floor(start / 0.1)) %>%
  summarize(mean_diversity = mean(diversity, na.rm = TRUE),
            start = first(start)) %>%
  ungroup()
ggplot(IS, aes(x = start, y = mean_diversity, group = group)) +
  geom_line(aes(color = group), size = 0.5) +
  labs(x = paste0("chromosome 12 position (Mb)"), 
       y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_color_manual(values=c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(marker~., scales='free')+
  scale_x_continuous(breaks = seq(34, 50, by = 1))+
  theme_bw()
ggsave(paste0(outPlot,"03.08.diversity_along-ImprovementSweep.pdf"),  width = 21, height = 20, units = "cm")  

  
#############################################
## D plot diversity for selective sweeps
#############################################
#sweeps<- data.table::fread("/home/IPS2/vibanez/Desktop/Q-lab/2024/03_dos-imp_sweeps/putative_dos-imp_sweeps.csv", sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
dir<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/ad_domestication-improvement-sweeps/"
sweeps<- data.table::fread(paste0(dir,"putative_dos-imp_sweeps.csv"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
toName<-c("chr", "start", "end",'1-WILD','2-PIM','3-SLC','4-SLL')
input<-list.files(pattern = "general_20acc_merged-pi.bed", full.names = FALSE)

diversity_sweep <- c()
noDiversity <- c()
diversity_wo_sweep <- c()
# Loop over input files
for(i in input) {
  df <- data.table::fread(i, sep = '\t', data.table = FALSE, fill = TRUE, na.string = "NA", nThread = 20)
  df$V8 <- NULL
  colnames(df)<-toName
  df$type <- unlist(strsplit(i, "_", fixed = TRUE))[1]
  # Subset by sweep
  sweep_results <- c()
  for(r in  1:nrow(sweeps)) {
    tmpSweep <- sweeps[r, ]
    tmp <- df[df$chr == tmpSweep[1, 2] & df$start >= tmpSweep[1, 3] & df$end <= tmpSweep[1, 4], ]
    result <- gsub("[0-9]", "", tmpSweep[1, 1])
    if (nrow(tmp) == 0) {
      noDiversity <- rbind.data.frame(noDiversity,tmpSweep[1, c(1, 2, 3, 4)])
    } else if (result == "DS") {
      tmp$sweep <- tmpSweep[1, 1]
      tmp$index <- tmp$`2-PIM` / tmp$`3-SLC`
    } else if (result == "IS") {
      tmp$sweep <- tmpSweep[1, 1]
      tmp$index <- tmp$`3-SLC` / tmp$`4-SLL`
    }
    sweep_results <- rbind.data.frame(sweep_results,tmp)
  }
  diversity_sweep<- rbind.data.frame(diversity_sweep,sweep_results)
  tmpDF <- anti_join(df, tmp, by = c("chr", "start", "end", "type"))
  diversity_wo_sweep <- rbind.data.frame(diversity_wo_sweep,tmpDF)
}

head(noDiversity)
head(diversity_sweep)
head(diversity_wo_sweep)

data.table::fwrite(diversity_sweep, file=paste0(outDir,"diversity_sweeep.tsv"), 
                   quote=F,row.names=F,col.names = T,sep="\t")

data.table::fwrite(diversity_wo_sweep, file=paste0(outDir,"diversity_wo_sweeep.tsv"), 
                   quote=F,row.names=F,col.names = T,sep="\t")

diversity_wo_sweep$sweep<-'nonSweep'
diversity_sweep$sweep<-gsub("[0-9]", "",diversity_sweep$sweep)
allDiversity<-rbind.data.frame(diversity_sweep[,c(infoCols,groupNames,'sweep','type')],
                               diversity_wo_sweep[,c(infoCols,groupNames,'sweep','type')])

allDiversity$chr<-gsub('SL2.50','',allDiversity$chr)
allDiversity <- reshape2::melt(allDiversity, id.vars = c("chr", "start", "end", "type", 'sweep'),
                               variable.name = "group",
                              value.name = "diversity")
allDiversity$diversity<-allDiversity$diversity*1000
quantiles <- quantile(allDiversity$diversity, probs = seq(0, 1, 0.01), na.rm = TRUE)
print(quantiles[91])
allDiversity$diversity[allDiversity$diversity >= 18] <- NA # missing for Q90= quantiles[91] | 90% |22.5907 
allDiversity <- na.omit(allDiversity)
head(allDiversity)
# only to check the index
{
summ<-plyr::ddply(allDiversity, c("type","sweep", "group"), summarise,
                  mean = mean(diversity, na.rm=TRUE)) #max = max(diversity, na.rm=TRUE),sd= sd(diversity, na.rm=TRUE),sem = sd(diversity,na.rm=TRUE)/sqrt(length(diversity))
tmpSumm<-reshape2::dcast(summ, type+sweep ~ group, value.var = "mean")
tmpSumm$WILD_PIM<-tmpSumm$`1-WILD`/tmpSumm$`2-PIM`
tmpSumm$PIM_SLC<-tmpSumm$`2-PIM`/tmpSumm$`3-SLC`
tmpSumm$SLC_SLL<-tmpSumm$`3-SLC`/tmpSumm$`4-SLL`
data.table::fwrite(tmpSumm, file=paste0(outDir,"diversity_index.tsv"), 
                   quote=F,row.names=F,col.names = T,sep="\t")
}
#allDiversity$samples_sweep <- paste(allDiversity$samples, allDiversity$type, sep = "_")

# plot
ggplot(allDiversity, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
  #geom_violin( na.rm = T) +
  geom_boxplot(width=0.5, na.rm = T,
               alpha=0.1,lwd = 0.3,outlier.shape = NA)+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) +
  labs(x = "",y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_fill_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  scale_color_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(sweep~type, scales = 'free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_line(color = "lightgray", linewidth = 0.2, linetype = "dashed"))+
  theme_bw()
ggsave(paste0(outPlot,"03.09.diversity_during_DOS-IMP.pdf"),  width = 20, height = 20, units = "cm")

tmp<-subset(allDiversity, group!="1-WILD")
ggplot(tmp, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
  #geom_violin( na.rm = T) +
  geom_boxplot(width=0.5, na.rm = T,
               alpha=0.1,lwd = 0.3,outlier.shape = NA)+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) +
  labs(x = "",y = expression(paste("(epi)", pi, " (10"^"-3"~")")))+
  scale_fill_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  scale_color_manual(values = c("#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(sweep~type, scales = 'free')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_line(color = "lightgray", linewidth = 0.2, linetype = "dashed"))+
  coord_cartesian(ylim = c(0, 12.5))+
  theme_bw()
ggsave(paste0(outPlot,"03.10.diversity_during_DOS-IMP_woWild.pdf"),  width = 20, height = 20, units = "cm")
}