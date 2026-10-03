######### load heritability
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(eulerr)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})
outDir<-"/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results/"
outPlot<-"/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/plots/"

## get the GIF of all sig DMRs:
allGIF<-fread("/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/af_get-nonSig-GIF/00.0_results-adjuted-lambda.tsv", 
              sep = '\t', header = T, fill = TRUE, na.strings = "",nThread = 40)
colnames(allGIF)<-c("DMR","dmr",'nAssociatedSNPs','lambda','sigFuncionLambda')
# h2 for all DMRs (sig+nonSig)
h2 <- fread(paste0(outDir,"00.1_merged-H2_allDMRs.tsv"), sep = '\t', header = T, 
            fill = TRUE, na.strings = "",nThread = 40)
setDT(h2)
# beta
betaVal<-fread(paste0(outDir,"00.1_merged-beta-values_sig-DMRs.tsv"), sep = '\t', header = T, 
               fill = TRUE, nThread = 40)
betaVal[, DMR := tstrsplit(DMRid, "_", fixed = TRUE)[2]]
betaVal$absValues<-abs(betaVal$beta)

## make the beta values
ggplot(betaVal, aes(x=DMR, y=absValues, fill=DMR, color=DMR)) +
  geom_violin(trim=FALSE, alpha=0.6) +  # Create the violin plot
  geom_boxplot(width=0.1, alpha=0.2, outlier.shape = NA, color="black", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#820a86", "#ffca7b")) +
  scale_color_manual(values=c("#820A86", "#ffca7b")) +
  labs(x = "DMR",y = "Absolute Beta Values") +
  theme_minimal() 
ggsave(paste0(outPlot,"00.0_beta_absValues.pdf"), width = 16, height = 12, units = "cm")

ggplot(betaVal, aes(x=DMR, y=beta, fill=DMR, color=DMR)) +
  geom_violin(trim=FALSE, alpha=0.6) +  # Create the violin plot
  geom_boxplot(width=0.1, alpha=0.2, outlier.shape = NA, color="black", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#820a86", "#ffca7b")) +
  scale_color_manual(values=c("#820A86", "#ffca7b"))+
  theme_minimal()
ggsave(paste0(outPlot,"00.0_beta_values.pdf"), width = 30, height = 20, units = "cm")

#setdiff(h2$DMR, allGIF$dmrPos)
# Perform the merge
#h2_merged <- merge(h2, maf, by = "DMR", all.x = TRUE)
h2_merged <- merge(h2, allGIF, by = "DMR", all.x = TRUE)
h2_merged <- h2_merged %>%
  mutate(associationType = case_when(
    type == "nonSig" & is.na(sigFuncionLambda) ~ "NonSig",
    type == "nonSig-GIF" & sigFuncionLambda == "linked_w_GI" ~ "NonSig_by_GIF",
    TRUE ~ 'sig_wo_GIF'
  ))

h2_merged[, c("chr", "DMRtype", "position") := tstrsplit(DMR, "_")]

h2_merged <- h2_merged %>%
  mutate(H2= sigmaGen / (sigmaGen + sigmaError) )#,
    #MAF_bin = cut(MAF, breaks = seq(0, 0.5, by = 0.05), include.lowest = TRUE))

h2_merged$associationType <- as.factor(h2_merged$associationType)
h2_merged$DMRtype <- as.factor(h2_merged$DMRtype)

ggplot(h2_merged, aes(x =H2, fill=associationType, color=associationType)) + 
  geom_density(alpha=0.1) +
  labs(title='',x="H2", y = "Density")+
  scale_fill_manual(values = c("#5ba48db2","#ccccccb2","darkred"))+
  scale_color_manual(values = c("#5ba48db2","#ccccccb2","darkred"))+
  xlim(0,1) +
  #geom_vline(data = mean_H2, aes(xintercept = meanH2, color = associationType), 
  #           linetype = "dashed", size = 0.8) +
  theme_minimal() +
  facet_grid(~DMRtype)
ggsave(paste0(outPlot,"00.1.H2_density.pdf"),  width = 16, height = 12, units = "cm")

#
plot_h2 <- h2_merged %>%
  filter(associationType %in% c("NonSig", "sig_wo_GIF")) %>%
  mutate(
    DMRtype = factor(DMRtype, levels = c("CG-DMR", "C-DMR")),
    link_status = case_when(
      associationType == "NonSig" ~ "Unlinked",
      associationType == "sig_wo_GIF"    ~ "Linked"
    ),
    link_status = factor(link_status, levels = c("Unlinked", "Linked"))
  )

ggplot(plot_h2, aes(x = DMRtype, y = H2, fill = DMRtype, color = DMRtype)) +
  geom_violin(alpha = 0.25, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "black") +
  scale_fill_manual(values = c("#ffca7b", "#820a86")) +
  scale_color_manual(values = c("#ffca7b", "#820A86")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "", x = "DMR type", y = "H2") +
  facet_wrap(~link_status) +
  theme_minimal()
ggsave(paste0(outPlot, "00.2.H2_violin_sig_only.pdf"),
       width = 16, height = 12, units = "cm")
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



