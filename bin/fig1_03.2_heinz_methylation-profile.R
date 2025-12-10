suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

setwd("/mnt/disk2/vibanez/10_data-analysis/Fig1/ab_global-methylation/heinz_methylomes/")
outPlot <- "/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
outDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"

# Load chromatic regions
chromatic_regions <- fread(paste0(outDir, "03.01_SL2.5_chromatin_regions.tsv"))
# Set bin size
window_size <- 2.5  # 2.5 Mb window
# Load all methylation files
input <- list.files(pattern = ".bed$", full.names = FALSE)
# Initialize
methylation_binned <- list()
for (i in input) {
  df <- fread(i, header = FALSE, col.names = c("chr", "start", "end", "meth"))
  acc <- strsplit(i, "_")[[1]][1]
  tissue <- strsplit(i, "_")[[1]][2]
  context <- gsub('.bed', '', strsplit(i, "_")[[1]][4])

  df$acc <- acc
  df$tissue <- tissue
  df$ctxt <- context
  
  # Convert start to Mb
  df$start_Mb <- df$start / 1e6
  
  # Bin
  df_binned <- df %>%
    mutate(window = floor(start_Mb / window_size)) %>%
    group_by(chr, ctxt, tissue, acc, window) %>%
    summarize(
      meanMeth = mean(meth, na.rm = TRUE),
      start_Mb = first(start_Mb),
      .groups = "drop"
    )
  
  methylation_binned[[i]] <- df_binned
}

# Combine all accessions
methylation_binned_df <- bind_rows(methylation_binned)
methylation_binned_df <- methylation_binned_df %>%
  mutate(
    start_bp = window * window_size * 1e6,     # window * 2.5Mb
    end_bp = start_bp + window_size * 1e6       # add window size
  )

# Add "SL2.50" prefix back
setnames(chromatic_regions, c("bin_start", "bin_end"), c("start", "end"))
# Set as data.table
setDT(methylation_binned_df)
setDT(chromatic_regions)

# Set keys
setkey(chromatic_regions, chr, start, end)
setkey(methylation_binned_df, chr, start_bp, end_bp)

# Overlap join
meth_with_region <- foverlaps(methylation_binned_df, chromatic_regions, nomatch = 0)

# Plot
ggplot() +
  # Background regions
  geom_rect(data = chromatic_regions,
            aes(xmin = start/1e6, xmax = end/1e6, ymin = -Inf, ymax = Inf, fill = region_type),
            alpha = 0.15) +
  # Smoothed methylation curves with linetype for tissue
  geom_smooth(data = meth_with_region,
              aes(x = start_bp/1e6, y = meanMeth, color = ctxt, linetype = tissue),
              method = "loess", span = 0.15, se = FALSE, size = 1) +
  facet_wrap(~chr, scales = "free_x", nrow = 3) +
  scale_color_manual(values=c("CG"= "#CCCC00","CHG"="#0033CC","CHH"= "#CC0066"))+
  scale_fill_manual(values = c("Euchromatic"="#66c2a5", "Heterochromatic"="#fc8d62", "Ambiguous"="white")) +
  scale_linetype_manual(values = c("leaf" = "solid", "fruit" = "dashed")) +  # Solid for leaf, dashed for fruit
  labs(x = "Genomic position (Mb)", 
       y = "Methylation level",
       fill = "Region type",
       color = "Context",
       linetype = "Tissue") +
  theme_bw() +
  ylim(0,1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10)
  )
ggsave(paste0(outPlot, "03.03.Heinz_methylation_genome-wide.pdf"), width = 51, height = 20, units = "cm")
##
chr1_meth <- meth_with_region %>%
  filter(chr == "SL2.50ch01")

chromatic_regions_ch1<-chromatic_regions%>%
  filter(chr == "SL2.50ch01")
  
ggplot() +
  # Background regions
  geom_rect(data = chromatic_regions_ch1,
            aes(xmin = start/1e6, xmax = end/1e6, ymin = -Inf, ymax = Inf, fill = region_type),
            alpha = 0.15) +
  # Smoothed methylation curves with linetype for tissue
  geom_smooth(data = chr1_meth,
              aes(x = start_bp/1e6, y = meanMeth, color = ctxt, linetype = tissue),
              method = "loess", span = 0.15, se = FALSE, size = 1) +
  #facet_wrap(~chr, scales = "free_x", nrow = 3) +
  scale_color_manual(values=c("CG"= "#CCCC00","CHG"="#0033CC","CHH"= "#CC0066"))+
  scale_fill_manual(values = c("Euchromatic"="#66c2a5", "Heterochromatic"="#fc8d62", "Ambiguous"="white")) +
  
  scale_linetype_manual(values = c("leaf" = "solid", "fruit" = "dashed")) +  # Solid for leaf, dashed for fruit

  labs(x = "Genomic position (Mb)", 
       y = "Methylation level",
       fill = "Region type",
       color = "Context",
       linetype = "Tissue") +
  
  theme_bw() +
  ylim(0,1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10)
  )
ggsave(paste0(outPlot, "03.04.Heinz_methylation_genome-wide_chr01.pdf"), width = 51, height = 20, units = "cm")

##
summary <- meth_with_region %>%
  group_by(region_type,tissue,ctxt, chr ) %>%
  summarise(mean_methylation = mean(meanMeth, na.rm = TRUE),
            .groups = "drop")
fwrite(summary, file = paste0(outDir, "03.02_Heinz_mean-meth_chromatin_regions.tsv"), sep = "\t", quote = FALSE)

summary$label<-paste0(summary$region_type, "_", summary$tissue)
ggplot(summary, aes(x = label, y = mean_methylation, fill = ctxt)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.75), alpha=0.5) +
  geom_jitter(aes(color = ctxt), 
              shape = 21, size = 1.5, 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75),
              alpha = 0.7) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = sprintf("%.2f", ..y..)),  # Format to 2 decimals
    vjust = -0.7,   # Move a little above the box
    position = position_dodge(width = 0.75),
    size = 3.5,     # Text size
    color = "black"
  ) +
  scale_fill_manual(values=c("CG"= "#CCCC00","CHG"="#0033CC","CHH"= "#CC0066"))+
  scale_color_manual(values=c("CG"= "#CCCC00","CHG"="#0033CC","CHH"= "#CC0066"))+
  labs(x = "Region type", y = "Mean methylation", fill = "Context") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  ylim(0, 1.1)  # Allow space for the labels above 1.0

ggsave(paste0(outPlot, "03.05.Heinz_methylation_chromatic_states.pdf"), width = 51, height = 20, units = "cm")
# A tibble: 18 × 4
# region_type     tissue ctxt  mean_methylation
# <chr>           <chr>  <chr>            <dbl>
# 1 Ambiguous       fruit  CG              0.828 
# 2 Ambiguous       fruit  CHG             0.673 
# 3 Ambiguous       fruit  CHH             0.165 
# 4 Ambiguous       leaf   CG              0.850 
# 5 Ambiguous       leaf   CHG             0.613 
# 6 Ambiguous       leaf   CHH             0.0831
# 7 Euchromatic     fruit  CG              0.573 
# 8 Euchromatic     fruit  CHG             0.170 
# 9 Euchromatic     fruit  CHH             0.102 
# 10 Euchromatic     leaf   CG              0.677 
# 11 Euchromatic     leaf   CHG             0.158   
# 12 Euchromatic     leaf   CHH             0.0549
# 13 Heterochromatic fruit  CG              0.813 
# 14 Heterochromatic fruit  CHG             0.652 
# 15 Heterochromatic fruit  CHH             0.163 
# 16 Heterochromatic leaf   CG              0.844 
# 17 Heterochromatic leaf   CHG             0.603 
# 18 Heterochromatic leaf   CHH             0.0857

