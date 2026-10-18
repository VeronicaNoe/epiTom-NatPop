suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(ggdendro) 
})

inDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/ba_formated-data/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/bb_sample-clustering/"
#########################
treeDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/aa_distances/"
treeOrder<- data.table::fread(paste0(treeDir,'general_SNPs.mdist.gz'), 
                              sep = '\t', data.table = F, 
                              fill = TRUE, na.string=c("NA"), nThread = 20) 
head(treeOrder)
colNames<-data.table::fread(paste0(treeDir, "general_SNPs.mdist.id"), 
                            sep = '\t', header = FALSE,
                            data.table = F, fill = TRUE, 
                            na.string=c("NA"), nThread = 20)[,2] 

# rename samples to match
colNames<-gsub('^0_','',colNames)
colNames<-gsub('S_','S-',colNames)
colNames<-gsub('_1','',colNames)
colNames<-gsub('BGV006775','S-pimBGV006775',colNames)
colnames(treeOrder)<-colNames
rownames(treeOrder)<-colNames
head(treeOrder)
# load sample names and groups
sampleTab<-data.table::fread("/mnt/disk2/vibanez/10_data-analysis/Fig1/00_sampleTab.tsv", 
                             sep = '\t', data.table = FALSE, fill = TRUE, 
                             na.string = "NA", nThread = 20)
sampleTab<-subset(sampleTab, Organ =="leaf")
head(sampleTab)

snps_clean<- data.table::fread(paste0(inDir,"SNP_leaf-metabolite_snp-wo-na-all-dmrs.tsv"), 
                                 sep = '\t', data.table = T, 
                               fill = TRUE, na.string=c("NA"), nThread = 20)
sampleNames <- snps_clean$Samples
subTree<-treeOrder[sampleNames,]
#
perform_mdist_clustering <- function(dist_matrix, k_clusters = 5, out_name = "genetic-distance-clustering", outDir) {
  # Ensure it's a matrix
  dist_matrix <- as.matrix(dist_matrix)
  # Convert distance matrix to 'dist' object
  dist_obj <- as.dist(dist_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_obj, method = "ward.D")
  clusters <- cutree(hc, k = k_clusters)
  
  # Save cluster table
  cluster_table <- data.frame(Sample = names(clusters), Cluster = clusters)
  outFile <- paste0(outDir, out_name, "_cluster.tsv")
  data.table::fwrite(cluster_table, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Dendrogram plot
  dendr_dendrogram <- as.dendrogram(hc)
  dendr_data <- dendro_data(dendr_dendrogram)
  
  # Add cluster info to labels
  dendr_data$labels$cluster <- as.factor(cluster_table$Cluster[match(dendr_data$labels$label, cluster_table$Sample)])
  
  # Plot
  plot <- ggplot() +
    geom_segment(data = dendr_data$segments, aes(x = x, y = y, xend = xend, yend = yend), size = 0.5) +
    geom_point(data = dendr_data$labels, aes(x = x, y = y, color = cluster), size = 2) +
    geom_text(data = dendr_data$labels, aes(x = x, y = y, label = label, color = cluster), hjust = 0, size = 2) +
    scale_color_manual(values = c("#FF7b89", "#8a5082", "#6f5f90", "#758eb7", "#a5cad2")) +
    theme_minimal() +
    coord_flip() +
    scale_y_reverse(expand = c(0.2, 0)) +
    ggtitle(paste("Dendrogram based on", out_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save plot
  ggsave(filename = paste0(outDir, out_name, "_dendrogram.pdf"), plot = plot, width = 20, height = 20, units = "cm")
}

perform_mdist_clustering(subTree, k_clusters = 5, out_name = "genetic-distance-clustering", outDir = outDir)

perform_mdist_clustering_by_group <- function(dist_matrix, sample_info, out_name = "genetic-distance-clustering", outDir) {
  # Ensure it's a matrix
  dist_matrix <- as.matrix(dist_matrix)
  
  # Convert distance matrix to 'dist' object
  dist_obj <- as.dist(dist_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_obj, method="ward.D")
  
  # Prepare group info
  sample_info <- sample_info[match(rownames(dist_matrix), sample_info$AccName), ]
  
  # Dendrogram plot
  dendr_dendrogram <- as.dendrogram(hc)
  dendr_data <- dendro_data(dendr_dendrogram)
  
  # Add Group info to labels
  dendr_data$labels$Group <- sample_info$Group[match(dendr_data$labels$label, sample_info$AccName)]
  
  # Pick colors (auto or custom)
  group_colors <- RColorBrewer::brewer.pal(length(unique(dendr_data$labels$Group)), "Set1")
  names(group_colors) <- unique(dendr_data$labels$Group)
  
  # Plot
  plot <- ggplot() +
    geom_segment(data = dendr_data$segments, aes(x = x, y = y, xend = xend, yend = yend), size = 0.5) +
    geom_point(data = dendr_data$labels, aes(x = x, y = y, color = Group), size = 2) +
    geom_text(data = dendr_data$labels, aes(x = x, y = y, label = label, color = Group), hjust = 0, size = 2) +
    scale_color_manual(values = group_colors) +
    theme_minimal() +
    coord_flip() +
    scale_y_reverse(expand = c(0.2, 0)) +
    ggtitle(paste("Dendrogram colored by Group (", out_name, ")", sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save plot
  ggsave(filename = paste0(outDir, out_name, "_dendrogram_by_group.pdf"), plot = plot, width = 20, height = 20, units = "cm")
}

perform_mdist_clustering_by_group(subTree, sampleTab, out_name = "genetic-distance-clustering", outDir = outDir)

