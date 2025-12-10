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

############
annotationDir<-"/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/"
genes<-data.table::fread(paste0(annotationDir,"allGenes.bed"), 
                        sep = '\t', data.table = T, header = F,
                        fill = TRUE, na.string=c("NA"), nThread = 20) 
infoCol<-c('chr','start','end','strand')
colnames(genes)<-c(infoCol,'name')
genes<-genes%>%
  filter(chr!="SL2.50ch00")%>% 
  mutate(type = 'gene')

TEs<-data.table::fread(paste0(annotationDir,"TE-wo-gene.anno"), 
                         sep = '\t', data.table = T, header = F,
                         fill = TRUE, na.string=c("NA"), nThread = 20) 
colnames(TEs)<-c(infoCol,'name')
TEs <- TEs %>% select(chr,start,end,strand, name)%>%
  filter(chr!="SL2.50ch00")%>%
  mutate(type = 'TE') 

####
chrom_sizes<-data.table::fread("~/bin/chr.size.bed", 
                         sep = '\t', data.table = T, header = F,
                         fill = TRUE, na.string=c("NA"), nThread = 20) 
colnames(chrom_sizes)<-c("chr","size")

bin_size <- 100000  # 100 kb
bins <- chrom_sizes %>%
  filter(chr!="SL2.50ch00") %>%
  rowwise() %>%
  do({
    chr = .$chr
    size = .$size
    bin_starts = seq(0, size, by = bin_size)
    data.frame(
      chr = chr,
      bin_start = bin_starts,
      bin_end = pmin(bin_starts + bin_size, size)
    )
  }) %>%
  ungroup()

# Function to compute overlap length
compute_overlap_length <- function(bin_chr, bin_start, bin_end, features_df) {
  features_chr <- features_df %>% filter(chr == bin_chr)
  
  # Find overlapping features
  overlaps <- features_chr %>%
    filter(start < bin_end & end > bin_start)
  
  if (nrow(overlaps) == 0) {
    return(0)
  }
  
  # Calculate overlap per feature
  overlap_lengths <- pmin(overlaps$end, bin_end) - pmax(overlaps$start, bin_start)
  overlap_lengths <- ifelse(overlap_lengths < 0, 0, overlap_lengths)  # just in case
  
  return(sum(overlap_lengths))
}
# Gene coverage
bins$gene_coverage <- mapply(
  compute_overlap_length,
  bin_chr = bins$chr,
  bin_start = bins$bin_start,
  bin_end = bins$bin_end,
  MoreArgs = list(features_df = genes)
)

# TE coverage
bins$TE_coverage <- mapply(
  compute_overlap_length,
  bin_chr = bins$chr,
  bin_start = bins$bin_start,
  bin_end = bins$bin_end,
  MoreArgs = list(features_df = TEs)
)


bins_plot <- bins %>%
  mutate(
    gene_fraction = gene_coverage / bin_size,
    TE_fraction = TE_coverage / bin_size,
    bin_midpoint = (bin_start + bin_end) / 2 / 1e6  # Midpoint in Mb
  )
# Colors
custom_colors <- c("blue", "#33a7cc", "yellow", "orange", "red")

# Plot gene coverage
plot_gene <- ggplot(bins_plot, aes(x = bin_midpoint, y = chr, fill = gene_fraction)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_colors, limits = c(0, 1)) +
  labs(title = "Gene Coverage (fraction of 100kb bin)",
       x = "Genomic position (Mb)",
       y = "Chromosome",
       fill = "Gene fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot TE coverage
plot_TE <- ggplot(bins_plot, aes(x = bin_midpoint, y = chr, fill = TE_fraction)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_colors, limits = c(0, 1)) +
  labs(title = "TE Coverage (fraction of 100kb bin)",
       x = "Genomic position (Mb)",
       y = "Chromosome",
       fill = "TE fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(paste0(outPlot,"03.01_genome_coverage_heatmap_gene_TE.pdf"), width = 12, height = 8)
print(plot_gene)
print(plot_TE)
dev.off()
##
bins_plot <- bins_plot %>%
  mutate(
    region_type = case_when(
      gene_fraction > TE_fraction ~ "Euchromatic",
      gene_fraction < TE_fraction ~ "Heterochromatic",
      TRUE ~ "Ambiguous"  # In case they are exactly equal (rare)
    )
  )

plot_regions <- ggplot(bins_plot, aes(x = bin_midpoint, y = chr, fill = region_type)) +
  geom_tile() +
  scale_fill_manual(values = c(
    "Euchromatic" = "#66c2a5",
    "Heterochromatic" = "#fc8d62",
    "Ambiguous" = "grey80"
  )) +
  labs(title = "Chromatin Type by Bin",
       x = "Genomic position (Mb)",
       y = "Chromosome",
       fill = "Region type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(paste0(outPlot,"03.02_chromatin_regions_by_bin.pdf"), width = 12, height = 8)
print(plot_regions)
dev.off()

fwrite(bins_plot, file = paste0(outDir, "03.01_SL2.5_chromatin_regions.tsv"), sep = "\t", quote = FALSE)

