suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggdendro)
  #library(ape)
  library(pheatmap)
  library(dendextend)
  library(dplyr)
  library(data.table)
  library(vegan)
  
})
#setwd("/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/data")
#outPlot<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/plots/"
#outDir<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/results/"
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig1/aa_distances")
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"
############################################################
##### dendrogram using 1-IBS distance
#################################################################
metaFile<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/00_sampleTab_remove.tsv"
sampleTab<-data.table::fread(metaFile, 
                             sep = '\t', data.table = FALSE, fill = TRUE, 
                             na.string=c("NA"), nThread = 16) 
sampleTab <- sampleTab %>%
  filter(Organ=="leaf")

input<-list.files(pattern = ".mdist.gz", full.names = FALSE)
for( i in input){
  sample<-unlist(strsplit(i,".", fixed=T))[1]
  mm<-unlist(strsplit(sample,"_", fixed=T))[2]
  dt<- data.table::fread(paste0(sample,".mdist.gz"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
  
  if(mm =="SNPs"){
    colNames<-data.table::fread(paste0(sample,".mdist.id"), sep = '\t', header = FALSE,
                                data.table = FALSE, fill = TRUE, 
                                na.string=c("NA"), nThread = 20)[,2] 
    colNames<-gsub("^0_","",colNames)
    colnames(dt)<-colNames
    rownames(dt)<-colNames
  }else {
    colNames<-data.table::fread(paste0(sample,".mdist.id"), sep = '\t', header = FALSE,
                                data.table = FALSE, fill = TRUE, 
                                na.string=c("NA"), nThread = 20)[,1] 
    colNames<-gsub("^0_","",colNames)
    colNames<-gsub("LA2838A","S-lycLA2838A",colNames)
    colnames(dt)<-colNames
    rownames(dt)<-colNames
  }
  dendr <- hclust(as.dist(dt))
  dendr_dendrogram <- as.dendrogram(dendr)
  # Convert dendrogram to a format compatible with ggplot2
  toPlot <- dendro_data(dendr_dendrogram)
  if(mm =="SNPs"){
    sampleTab <- sampleTab %>%
      mutate(SNPnames = ifelse(SNPnames == "", AccName, SNPnames))
    toKeep <- intersect(toPlot$labels$label, sampleTab$SNPnames)
    rownames(sampleTab) <- sampleTab$SNPnames
    group <- sampleTab[toKeep, "Group"]
    toPlot$labels$group <- group
  }else{
    toKeep <- intersect(toPlot$labels$label, sampleTab$BiseqName)
    rownames(sampleTab) <- sampleTab$BiseqName
    group <- sampleTab[toKeep, "Group"]
    toPlot$labels$group <- group
  }
  # Define color map
  color_map <- c("#000000ff","#00a681ff","#ed83b5ff","#00a4deff")
  names(color_map) <- c("WILD","PIM","SLC","SLL")
  toPlot$labels$group <- factor(toPlot$labels$group, levels = c("WILD", "PIM", "SLC", "SLL"))
  # Plot the dendrogram
  ggplot() +
    geom_segment(data = toPlot$segments,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(data = toPlot$labels,
               aes(x = x, y = y, color = group),
               size = 2) +
    geom_text(data = toPlot$labels,
              aes(x = x, y = y, label = label, hjust = 0, colour = group), hjust = 0, size = 1) +
    scale_color_manual(values = color_map) +
    theme_minimal() +
    coord_flip() +
    scale_y_reverse(expand = c(0.2, 0))
  # Save the plot
  ggsave(paste0(outPlot,"01.01_dendrogram_",sample,"_1-ibs.pdf"), width = 20, height = 20, units = "cm")
  
  #get order
  orderLeaf<-as.data.frame(toPlot$labels)
  write.table(orderLeaf, paste0(outDir,"01.01_data_",sample,"_order.csv"),quote=F,row.names=F,sep="\t")
}

############################################################
##### correlation between distance matrices
#################################################################
# Read sample table
sampleTab <- fread(metaFile, sep = '\t', data.table = FALSE, fill = TRUE, na.string = c("NA"), nThread = 16) %>%
  filter(Organ == "leaf") %>%
  mutate(SNPnames = ifelse(SNPnames == "", AccName, SNPnames))

rownames(sampleTab) <- sampleTab$SNPnames

# Read SNP distance matrix
snpDis <- fread("general_SNPs.mdist.gz", sep = '\t', data.table = FALSE, fill = TRUE, na.string = c("NA"), nThread = 20)
snpNames <- fread("general_SNPs.mdist.id", sep = '\t', header = FALSE, data.table = FALSE, fill = TRUE, na.string = c("NA"), nThread = 20)[,2]
snpNames <- gsub("^0_", "", snpNames)
snpSamples <- sampleTab[snpNames,]
colnames(snpDis) <- snpSamples$AccName
rownames(snpDis) <- snpSamples$AccName

# Initialize result tables
correlation_results <- data.frame(DMR_type = character(), Pearson_cor = numeric(), Pearson_pval = numeric())
mantel_results <- data.frame(DMR_type = character(), Mantel_r = numeric(), Mantel_pval = numeric())
cophenetic_results <- data.frame(DMR_type = character(), Cophenetic_cor = numeric())
entanglement_results <- data.frame(DMR_type = character(), Entanglement = numeric())

# DMR types to analyze
dmr_types <- c("CG-DMR", "C-DMR")

for (i in dmr_types) {
  cat(paste0("\nProcessing: ", i, "\n"))
  
  # Read DMR distance matrix
  dmrDis <- fread(paste0("general_", i, ".mdist.gz"), sep = '\t', data.table = FALSE, fill = TRUE, na.string = c("NA"), nThread = 20)
  dmrNames <- fread(paste0("general_", i, ".mdist.id"), sep = '\t', header = FALSE, data.table = FALSE, fill = TRUE, na.string = c("NA"), nThread = 20)[,1]
  dmrNames <- gsub("LA2838A", "S-lycLA2838A", dmrNames)
  rownames(sampleTab) <- sampleTab$BiseqName
  dmrSamples <- sampleTab[dmrNames,]
  colnames(dmrDis) <- dmrSamples$AccName
  rownames(dmrDis) <- dmrSamples$AccName
  
  # Reorder matrices
  common_names <- intersect(rownames(dmrDis), rownames(snpDis))
  dmrDis <- dmrDis[common_names, common_names]
  snpDis <- snpDis[common_names, common_names]
  
  # Vectors for correlation
  dmr_vector <- as.vector(dmrDis[lower.tri(dmrDis)])
  snp_vector <- as.vector(snpDis[lower.tri(snpDis)])
  
  # Pearson correlation
  pearson_result <- cor.test(dmr_vector, snp_vector, method = "pearson")
  correlation_results <- rbind(correlation_results,
                               data.frame(DMR_type = i,
                                          Pearson_cor = pearson_result$estimate,
                                          Pearson_pval = pearson_result$p.value))
  
  # Mantel test
  mantel_result <- mantel(as.dist(snpDis), as.dist(dmrDis), method = "pearson", permutations = 9999)
  mantel_results <- rbind(mantel_results,
                          data.frame(DMR_type = i,
                                     Mantel_r = mantel_result$statistic,
                                     Mantel_pval = mantel_result$signif))
  
  # Dendrograms and cophenetic correlation
  dmrDend <- as.dendrogram(hclust(as.dist(dmrDis)))
  snpDend <- as.dendrogram(hclust(as.dist(snpDis)))
  cophenetic_correlation <- cor_cophenetic(dmrDend, snpDend)
  cophenetic_results <- rbind(cophenetic_results,
                              data.frame(DMR_type = i,
                                         Cophenetic_cor = cophenetic_correlation))
  
  # Better untangling
  untangled <- untangle(dmrDend, snpDend, method = "step2side")
  entanglement_value <- entanglement(untangled)
  entanglement_results <- rbind(entanglement_results,
                                data.frame(DMR_type = i,
                                           Entanglement = entanglement_value))
  
  # First, your manual color mapping:
  group_colors <- c("WILD" = "#000000ff", 
                    "PIM" = "#00a681ff", 
                    "SLC" = "#ed83b5ff", 
                    "SLL" = "#00a4deff")
  
  # Match your samples to their groups
  common_samples <- intersect(labels(dmrDend), sampleTab$AccName)
  sampleTab_common <- sampleTab[match(common_samples, sampleTab$AccName), ]
  
  # Now assign the manual colors
  label_colors <- group_colors[sampleTab_common$Group]
  names(label_colors) <- sampleTab_common$AccName
  
  # Untangle the dendrograms
  dend_list <- untangle(dendlist(dmrDend, snpDend), method = "step2side")
  # Make sure label_colors are ordered correctly
  ordered_colors1 <- label_colors[labels(dend_list[[1]])]
  ordered_colors2 <- label_colors[labels(dend_list[[2]])]
  
  # Now color_labels with correct matching
  dend_list[[1]] <- color_labels(dend_list[[1]], col = ordered_colors1)
  dend_list[[2]] <- color_labels(dend_list[[2]], col = ordered_colors2)
  
  
  # Plot
  pdf(paste0(outPlot,"01.02_tanglegram_",i,".pdf"), width = 10, height = 8)
  tanglegram(dend_list[[1]], dend_list[[2]],
             highlight_distinct_edges = FALSE,
             common_subtrees_color_lines = FALSE,
             lwd = 1,
             color_lines = "grey",
             lab.cex = 0.6,
             margin_inner = 10,
             main = "DMRs vs SNPs Tanglegram")
  dev.off()
  
  # QQ plots
  pdf(paste0(outPlot, "01.03_qqplot_", i, "_SNPs.pdf"))
  par(mfrow = c(1,2))
  qqnorm(dmr_vector, main = paste0("QQ Plot of ", i, " Distances"))
  qqline(dmr_vector)
  qqnorm(snp_vector, main = "QQ Plot of SNP Distances")
  qqline(snp_vector)
  dev.off()
}

# Merge all results
overall_summary <- merge(correlation_results, mantel_results, by = "DMR_type") %>%
  merge(., cophenetic_results, by = "DMR_type") %>%
  merge(., entanglement_results, by = "DMR_type")

# Save result tables
write.table(overall_summary, paste0(outDir, "01.02_correlation_results.tsv"), sep = '\t', row.names = FALSE, quote = FALSE)

# Reshape for plotting
plot_data <- overall_summary %>%
  pivot_longer(cols = c(Pearson_cor, Mantel_r, Cophenetic_cor, Entanglement),
               names_to = "Metric", values_to = "Value")

# Plot
summary_plot <- ggplot(plot_data, aes(x = DMR_type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Summary of DMR vs SNP Comparisons", y = "Value", x = "DMR Type") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
pdf(paste0(outPlot,"01.04_summary_barplot_DMR_SNP_comparisons.pdf"), width = 8, height = 6)
print(summary_plot)
dev.off()
