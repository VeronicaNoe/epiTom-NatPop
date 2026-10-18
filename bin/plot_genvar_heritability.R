#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# =========================
# PATHS
# =========================
basedir <- "/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS/ag_outputs/"

asso_dirs <- c("sig", "nonSig")

# =========================
# HELPERS
# =========================
parse_file <- function(file) {
  stem <- sub("\\.QTL$", "", file)
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  
  data.table(
    metabolite = parts[1],
    mm         = parts[2],
    kinship    = parts[3],
    index      = parts[4],
    tag        = stem
  )
}

read_reml_h2 <- function(reml_file) {
  if (!file.exists(reml_file)) {
    warning("Missing REML file: ", reml_file)
    return(list(h2 = NA_real_, H2 = NA_real_))
  }
  
  x <- fread(reml_file, sep = "\t", fill = TRUE, header = FALSE)
  v <- suppressWarnings(as.numeric(x[[1]]))
  
  if (length(v) < 6) {
    warning("Unexpected REML format: ", reml_file)
    return(list(h2 = NA_real_, H2 = NA_real_))
  }
  
  sigmaGen   <- v[4]
  sigmaError <- v[5]
  h2_emmax   <- v[6]
  
  H2_manual <- sigmaGen / (sigmaGen + sigmaError)
  
  list(h2 = h2_emmax, H2 = H2_manual)
}

# =========================
# READ GENVAR RESULTS
# =========================
general_info <- list()

for (asso in asso_dirs) {
  resultType <- asso
  qtl_dir <- file.path(basedir, asso)
  files <- list.files(qtl_dir, pattern = "\\.QTL$", full.names = FALSE)
  for (file in files) {
    meta <- parse_file(file)
    # keep only DMR with genvar kinship
    if (!(meta$mm == "DMR" && meta$kinship == "genvar")) next
    qtl_file  <- file.path(qtl_dir, file)
    reml_file <- sub("\\.QTL$", ".reml", qtl_file)
    h2_vals <- read_reml_h2(reml_file)
    if (resultType == "sig") {
      qtl <- fread(qtl_file, sep = "\t", fill = TRUE, header = FALSE)
      nLoci <- nrow(qtl)
    } else {
      nLoci <- 0
    }
    general_info[[length(general_info) + 1]] <- data.table(
      metabolite = meta$metabolite,
      mm         = meta$mm,
      kinship    = meta$kinship,
      index      = meta$index,
      resultType = resultType,
      nLoci      = nLoci,
      h2         = h2_vals$h2,
      H2         = h2_vals$H2,
      tag        = meta$tag
    )
  }
}

general_info <- rbindlist(general_info, fill = TRUE)

# =========================
# SAVE TABLE
# =========================
fwrite(general_info, file.path(basedir, "genvar_heritability_general_info.tsv"),  sep = "\t",  quote = FALSE)

cat("\n=== GENVAR SUMMARY ===\n")
print(general_info %>%
    group_by(index) %>%
    summarise(
      n_traits = n(),
      mean_h2 = mean(h2, na.rm = TRUE),
      median_h2 = median(h2, na.rm = TRUE),
      mean_H2 = mean(H2, na.rm = TRUE),
      median_H2 = median(H2, na.rm = TRUE),
      .groups = "drop"))

# =========================
# PLOT 1: DENSITY H2 BY INDEX AND RESULT TYPE
# =========================
plot_df <- general_info %>%
  filter(is.finite(h2)) %>%
  mutate(index = factor(index, levels = c("aBN", "aIBS")),
         resultType = factor(resultType, levels = c("sig", "nonSig"))  )

summary_df <- plot_df %>%
  group_by(index) %>%
  summarise(mean_h2 = mean(h2, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(plot_df, aes(x = h2, color = index, fill = index)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  geom_vline(
    data = summary_df,
    aes(xintercept = mean_h2, color = index),
    linetype = "dashed",
    linewidth = 0.7  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),    limits = c(0, 1)  ) +
  scale_color_manual(values = c("aBN" = "#820a86", "aIBS" = "#999999")) +
  scale_fill_manual(values = c("aBN" = "#820a86", "aIBS" = "#999999")) +
  labs(x = "Heritability", y = "Density", title = "DMR metabolite GWAS heritability using SNP+TIP+SV genvar kinship"  ) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "white", color = NA),legend.title = element_blank()  )

ggsave(  file.path(basedir, "genvar_heritability_density_by_resultType.pdf"), p1,  width = 18,  height = 10,  units = "cm")
