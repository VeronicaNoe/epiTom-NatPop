########## get the MAF
From other FIGURE
{
# suppressPackageStartupMessages({
#   library(data.table)
#   library(dplyr)
#   library(eulerr)
#   library(dplyr)
#   library(ggplot2)
# })
# setwd("/home/IPS2/vibanez/Desktop/Q-lab/2023/SFS/data")
# outDir<-"/home/IPS2/vibanez/Desktop/Q-lab/2024/06_DMR-GWAS/results/"
# 
# # Read MAF data
# input<-list.files(pattern = "_SFS.tsv", full.names = FALSE)
# input<-grep('SNPs',input, value=T, invert = T)
# 
# maf_list <- lapply(input, function(file_path) {
#   out <- fread(file_path, sep = '\t', header = T, fill = TRUE, na.strings = "",nThread = 40)
#   return(out)
# })
# maf <- rbindlist(maf_list)
# colnames(maf)<-c('toRemove','toRemove','toRemove','toRemove','MAF','toRemove','toRemove',
#                 'DMR')
# maf<-maf[,!'toRemove']
# data.table::fwrite(maf, file= paste0(outDir,'00_merged-maf.tsv'), quote=F,
#                    row.names=F,col.names = T, sep="\t")
}
######### load heritability
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(eulerr)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/plots/"

## get the GIF of all sig DMRs:
allGIF<-fread(paste0(outDir,"00.0_results-adjuted-lambda.tsv"), sep = '\t', header = T, 
              fill = TRUE, na.strings = "",nThread = 40)
colnames(allGIF)<-c("DMR","dmr",'nAssociatedSNPs','lambda','sigFuncionLambda')
# h2 for all DMRs (sig+nonSig)
h2 <- fread(paste0(outDir,"00.1_merged-H2_allDMRs.tsv"), sep = '\t', header = T, 
            fill = TRUE, na.strings = "",nThread = 40)
setDT(h2)
# maf for all DMRs (common and nonCommon)
FROM OTHER FIGURE
maf <- fread(paste0(outDir,"00.0_merged-maf.tsv"), sep = '\t', header = T, 
            fill = TRUE, na.strings = "",nThread = 40) # all even those MAF<0.05
setDT(maf)
# beta
betaVal<-fread(paste0(outDir,"00.1_merged-beta-values_sig-DMRs.tsv"), sep = '\t', header = F, 
               fill = TRUE, nThread = 40)
betaVal[, DMR := tstrsplit(DMRpos, "_", fixed = TRUE)[2]]
betaVal$absValues<-abs(betaVal$beta)

## make the beta values
ggplot(betaVal, aes(x=DMR, y=absValues, fill=DMR, color=DMR)) +
  geom_violin(trim=FALSE, alpha=0.6) +  # Create the violin plot
  geom_boxplot(width=0.1, alpha=0.2, outlier.shape = NA, color="black", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#820a86", "#ffca7b")) +
  scale_color_manual(values=c("#820A86", "#ffca7b")) +
  labs(x = "DMR",y = "Absolute Beta Values") +
  theme_minimal() 
ggsave(paste0(outPlot,"00.0_beta_absValues.pdf"), width = 30, height = 20, units = "cm")

ggplot(betaVal, aes(x=DMR, y=beta, fill=DMR, color=DMR)) +
  geom_violin(trim=FALSE, alpha=0.6) +  # Create the violin plot
  geom_boxplot(width=0.1, alpha=0.2, outlier.shape = NA, color="black", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#820a86", "#ffca7b")) +
  scale_color_manual(values=c("#820A86", "#ffca7b"))+
  theme_minimal()
ggsave(paste0(outPlot,"00.0_beta_values.pdf"), width = 30, height = 20, units = "cm")

#setdiff(h2$DMR, allGIF$dmrPos)
# Perform the merge
h2_merged <- merge(h2, maf, by = "DMR", all.x = TRUE)
h2_merged <- merge(h2_merged, allGIF, by = "DMR", all.x = TRUE)
h2_merged<- h2_merged %>%
  mutate(associationType = case_when(
    type == "NonSig" & is.na(sigFuncionLambda) ~ "NonSig",
    type == "NonSig" & sigFuncionLambda == "linked_w_GI" ~ "NonSig_by_GIF",
    TRUE ~ 'sig_wo_GIF'
  ))

h2_merged[, c("chr", "DMRtype", "position") := tstrsplit(DMR, "_")]

h2_merged <- h2_merged %>%
  mutate(H2= sigmaGen / (sigmaGen + sigmaError),
    MAF_bin = cut(MAF, breaks = seq(0, 0.5, by = 0.05), include.lowest = TRUE))

h2_merged$associationType <- as.factor(h2_merged$associationType)
h2_merged$DMRtype <- as.factor(h2_merged$DMRtype)

# Create the density plot
{#H2-MAF
ggplot(h2_merged, aes(x = H2, y = MAF_bin, fill = associationType, 
                      color=associationType)) +
  geom_density_ridges(scale = 1, alpha = 0.8, alpha=0.4) +
  facet_grid(.~ DMRtype ) +
  scale_fill_manual(values = c("grey70","grey40","grey10"))+
  scale_color_manual(values = c("grey70","grey40","grey10"))+
  labs(title = "Density Plot of Heritability (H2) vs MAF Bin",
       x = "Heritability (H2)",
       y = "MAF Bin",
       fill = "Type") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outPlot,"00.1.H2-MAF_density.pdf"),  width = 40, height = 20, units = "cm")

ggplot(h2_merged, aes(x = MAF_bin, y = H2, color=associationType, fill=associationType)) +
  geom_boxplot(width=0.5, na.rm = T,alpha=0.6,lwd = 0.3,outlier.shape = NA)+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred",
               position = position_dodge(width = 0.5)) +
  #geom_violin(width = 0.7, alpha = 0.4, lwd = 0.3, position = position_dodge(width = 0.75)) +
  #stat_compare_means() +   # Add global p-value
  scale_fill_manual(values = c("grey70","grey40","grey10"))+
  scale_color_manual(values = c("grey70","grey40","grey10"))+
  labs(title = "Boxplot of h2 vs MAF Bins", x = "MAF Bins", y = "Heritability (h2)") +
  theme_bw() +
  facet_grid(associationType ~ DMRtype, scales = 'free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "bottom") +
  theme(panel.grid.major = element_line(color = "gray90", linetype = "dashed"))
ggsave(paste0(outPlot,"00.1.H2-MAF_boxplot.pdf"),  width = 40, height = 20, units = "cm")
}

mean_H2 <- h2_merged %>%
  group_by(associationType,  DMRtype) %>%
  summarize(meanH2 = mean(H2, na.rm = TRUE))
# associationType DMRtype meanH2
# <fct>           <fct>    <dbl>
# 1 NonSig          C-DMR    0.440
# 2 NonSig          CG-DMR   0.389
# 3 NonSig_by_GIF   C-DMR    0.557
# 4 NonSig_by_GIF   CG-DMR   0.538
# 5 sig_wo_GIF      C-DMR    0.629
# 6 sig_wo_GIF      CG-DMR   0.565
ggplot(h2_merged, aes(x =H2, fill=associationType, color=associationType)) + 
  geom_density(alpha=0.1) +
  labs(title='',x="H2", y = "Density")+
  scale_fill_manual(values = c("#5ba48db2","#ccccccb2","darkred"))+
  scale_color_manual(values = c("#5ba48db2","#ccccccb2","darkred"))+
  xlim(0,1) +
  geom_vline(data = mean_H2, aes(xintercept = meanH2, color = associationType), 
             linetype = "dashed", size = 0.8) +
  theme_minimal() +
  facet_grid(~DMRtype)
ggsave(paste0(outPlot,"00.1.H2_density.pdf"),  width = 40, height = 20, units = "cm")


### pieplot
summary_data <- h2_merged %>%
  filter(!is.na(DMRtype) & !is.na(associationType)) %>%
  group_by(DMRtype, associationType) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(DMRtype) %>%
  mutate(percentage = round(count / sum(count) * 100,1)) %>%
  ungroup()

# Function to create pie plots
create_pie_plot <- function(data, dmr) {
  ggplot(data, aes(x = "", y = count, fill = associationType)) +
    geom_bar(stat = "identity", width = 1) +  
    geom_text(aes(label = paste0(round(percentage, 2), "%")), position = position_stack(vjust = 0.5)) +  # Add labels inside slices
    coord_polar("y", start = 0) +  # Make the chart polar
    theme_void() +  # Remove unnecessary elements
    scale_fill_manual(values = c("#5ba48db2","#ccccccb2","darkred"))+
    labs(title = paste("Pie chart for", dmr), fill = "Association Type")
}
# Create and save the plots
unique_dmrs <- unique(summary_data$DMRtype)
pdf(paste0(outPlot,"00.2_DMRtype_associationType_pie.pdf"))
for (dmr in unique_dmrs) {
  data <- summary_data %>% filter(DMRtype == dmr)
  print(create_pie_plot(data, dmr))
}
dev.off()
# do a bar plot
ggplot(summary_data, aes(x=DMRtype, y=percentage, fill=associationType)) +
  geom_bar(stat="identity", position="stack") +  # Stacked bar plot
  scale_fill_manual(values = c("#5ba48db2", "#ccccccb2", "darkred")) +  # Color scale
  labs(x = "DMR Type", y = "Count", fill = "Association Type") +  # Labels
  theme_minimal()
ggsave(paste0(outPlot,"00.2_DMRtype_associationType_bar.pdf"), width = 30, height = 20, units = "cm")


