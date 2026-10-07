suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

in_dir <- "/mnt/disk2/vibanez/10_data-analysis/Fig1/ac_pca/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
markers <- c("CG-DMR", "C-DMR", "SNPs")

# Initialize list to store scree data
scree_list <- list()

for (m in markers) {
  # Load eigenvalues using scan (1 column numeric)
  eig <- scan(file.path(in_dir, paste0("general_", m, ".eigenval")))
  # Build scree data frame
  scree_df <- data.frame(
    Component = seq_along(eig),
    Variance = eig / sum(eig) * 100,
    Cumulative = cumsum(eig / sum(eig) * 100),
    Dataset = m
  )
  # Store in list
  scree_list[[m]] <- scree_df
}

# Combine all for plotting
scree_data <- do.call(rbind, scree_list)
#
threshold_df <- scree_data %>%
  group_by(Dataset) %>%
  filter(Cumulative >= 90) %>%
  slice_min(Component, n = 1) %>%
  ungroup()

#
# Run KS tests
ks_cg <- ks.test(scree_data$Cumulative[scree_data$Dataset == "SNPs"],
                 scree_data$Cumulative[scree_data$Dataset == "CG-DMR"])

ks_c <- ks.test(scree_data$Cumulative[scree_data$Dataset == "SNPs"],
                scree_data$Cumulative[scree_data$Dataset == "C-DMR"])

# Build summary table
ks_table <- data.frame(
  Comparison = c("SNPs vs CG-DMR", "SNPs vs C-DMR"),
  D = c(ks_cg$statistic, ks_c$statistic),
  p_value = c(ks_cg$p.value, ks_c$p.value)
)
ks_table
# Comparison    D    p_value
# 1 SNPs vs CG-DMR 0.45 0.03484457
# 2  SNPs vs C-DMR 0.35 0.17247627


# Plot
ggplot(scree_data, aes(x = Component, y = Cumulative, color = Dataset)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray40", size = 1) +
  geom_point(data = threshold_df, aes(x = Component, y = Cumulative), 
             shape = 21, fill = "white", size = 3, stroke = 1.2) +
  geom_text(data = threshold_df, aes(x = Component, y = Cumulative + 3, 
                                     label = paste0("PC", Component)), 
            color = "black", size = 3.5) +
  annotate("text", x = 3, y = 25, hjust = 0, size = 4, color = "black",
           label = paste0("KS Dstat and p-values:\n",
                          "SNPs_vs_CG-DMR: D = ", formatC(ks_cg$statistic, digits = 2, format = "f"), " p = ", formatC(ks_cg$p.value, digits = 3, format = "f"), "\n",
                          "SNPs_vs_C-DMR: D = ", formatC(ks_c$statistic, digits = 2, format = "f") ," p = ", formatC(ks_c$p.value, digits = 3, format = "f"))) +

  annotate("text", x = max(scree_data$Component) * 0.8, y = 91.5, 
           label = "90% Variance", color = "gray40", size = 4, hjust = 0) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b", "SNPs" = "#999999")) +
  labs(title = "Cumulative Explained Variance (PLINK PCA)",
       x = "Principal Component",
       y = "Cumulative Variance Explained (%)",
       color = "Marker Type") +
  ylim(0, 100)

ggsave(paste0(outPlot, "06.01_eigenvalues_cumulative.pdf"), width = 20, height = 20, units = "cm")

ggplot(scree_data, aes(x = Component, y = Variance, color = Dataset)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=c("C-DMR"="#820a86", "CG-DMR"="#ffca7b", "SNPs"="#999999"))+
  labs(title = "Explained Variance (PLINK PCA)",
       x = "Principal Component",
       y = "Variance Explained (%)",
       color = "Marker Type") +
  ylim(0, 50)
ggsave(paste0(outPlot,"06.01_eigenvalues_variance.pdf"), width = 20, height = 20, units = "cm")