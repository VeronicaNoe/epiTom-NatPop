#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# =========================
# PATHS
# =========================
workDir <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/am_pairwise-divergences"
outDir  <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/results/"
outPlot <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/plots/"

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
dir.create(outPlot, recursive = TRUE, showWarnings = FALSE)

setwd(workDir)

# =========================
# PARAMETERS
# =========================
DIFF_SCALE <- 1000

color_map <- c(
  "G0:G0" = "#00a4deff",
  "G0:G5" = "#a900de",
  "G5:G5" = "#de3a00"
)

# =========================
# HELPERS
# =========================
parse_sdiff_name <- function(file) {
  stem <- sub("\\.sdiff\\.summary$", "", basename(file))
  parts <- strsplit(stem, "\\.")[[1]]
  
  # SNP file: general.SNP-denovo.GQ20_DP10.sdiff.summary
  if (grepl("^SNP", parts[2])) {
    return(list(
      marker = "SNP",
      tech = "ONT",
      value_col = "SNPs-ONT",
      n_col = "nSNP-ONT"
    ))
  }
  
  # DMR files: general.BIseq.C-DMR-denovo.sdiff.summary
  tech <- parts[2]
  marker <- sub("-denovo.*$", "", parts[3])
  
  value_col <- paste(marker, tech, sep = "-")
  n_col <- paste0("n", value_col)
  
  list(
    marker = marker,
    tech = tech,
    value_col = value_col,
    n_col = n_col
  )
}

read_sdiff <- function(file) {
  info <- parse_sdiff_name(file)
  
  dt <- fread(file, sep = "\t", fill = TRUE, nThread = 20)
  
  if ("#IID1" %in% names(dt)) {
    setnames(dt, "#IID1", "IID1")
  }
  
  needed <- c("IID1", "IID2", "OBS_CT", "DIFF_CT")
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", file, ": ", paste(missing_cols, collapse = ", "))
  }
  
  dt <- dt[, ..needed]
  
  # Standardize pair orientation
  dt[, acc1 := fifelse(IID1 < IID2, IID1, IID2)]
  dt[, acc2 := fifelse(IID1 < IID2, IID2, IID1)]
  dt[, sampleComparison := paste(acc1, acc2, sep = ":")]
  
  dt[, DIFF_CT := DIFF_CT / DIFF_SCALE]
  
  setnames(dt, c("OBS_CT", "DIFF_CT"), c(info$n_col, info$value_col))
  
  dt[, .SD, .SDcols = c("sampleComparison", "acc1", "acc2", info$n_col, info$value_col)]
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

create_lm_plot <- function(data, x_var, y_var, color_by, output_file,
                           x_label = x_var, y_label = y_var,
                           outPlot, outDir, color_map) {
  
  df <- data %>%
    filter(
      !is.na(.data[[x_var]]),
      !is.na(.data[[y_var]]),
      !is.na(.data[[color_by]])
    ) %>%
    arrange(.data[[x_var]])
  
  if (nrow(df) < 3) {
    warning("Skipping plot; fewer than 3 rows: ", y_var, " vs ", x_var)
    return(NULL)
  }
  
  fit <- lm(as.formula(paste0("`", y_var, "` ~ 0 + `", x_var, "`")), data = df)
  
  newdata <- data.frame(
    x_tmp = seq(
      min(df[[x_var]], na.rm = TRUE),
      max(df[[x_var]], na.rm = TRUE),
      length.out = 200
    )
  )
  colnames(newdata) <- x_var
  
  pred <- as.data.frame(predict(fit, newdata = newdata, interval = "confidence"))
  predicted <- cbind(newdata, pred)
  colnames(predicted)[colnames(predicted) == "fit"] <- y_var
  
  fwrite(
    as.data.table(predicted),
    file.path(outDir, paste0(safe_name(y_var), "_vs_", safe_name(x_var), "_lm_predicted_values.tsv")),
    sep = "\t"
  )
  
  slope <- coef(fit)[1]
  r2 <- summary(fit)$r.squared
  pval <- summary(fit)$coefficients[1, 4]
  
  lab_txt <- paste0(
    "y = ", round(slope, 4), "x\n",
    "R² = ", round(r2, 3), "\n",
    "p = ", signif(pval, 3)
  )
  
  max_y <- max(df[[y_var]], na.rm = TRUE)
  min_x <- min(df[[x_var]], na.rm = TRUE)
  max_x <- max(df[[x_var]], na.rm = TRUE)
  
  p_main <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = .data[[color_by]]), alpha = 0.7, size = 2.5) +
    geom_ribbon(
      data = predicted,
      aes(x = .data[[x_var]], ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      fill = "grey70",
      alpha = 0.35
    ) +
    geom_line(
      data = predicted,
      aes(x = .data[[x_var]], y = .data[[y_var]]),
      inherit.aes = FALSE,
      color = "blue",
      linetype = "dashed",
      linewidth = 0.9
    ) +
    scale_color_manual(values = color_map, drop = FALSE) +
    annotate(
      "text",
      x = min_x + 0.05 * (max_x - min_x),
      y = max_y * 0.95,
      label = lab_txt,
      hjust = 0,
      vjust = 1,
      size = 4,
      fontface = "bold"
    ) +
    labs(x = x_label, y = y_label, color = "Comparison") +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  ggsave(
    file.path(outPlot, output_file),
    plot = p_main,
    width = 20,
    height = 20,
    units = "cm"
  )
  
  list(model = fit, predicted = predicted, plot = p_main)
}

# =========================
# 1) LOAD ALL SDIFF FILES
# =========================
input <- list.files(pattern = "\\.sdiff\\.summary$", full.names = FALSE)

processed_files <- lapply(input, read_sdiff)

merged_data <- Reduce(
  function(x, y) merge(x, y, by = c("sampleComparison", "acc1", "acc2"), all = TRUE),
  processed_files
)

# group label
merged_data[, group := paste0(
  sub("-.*$", "", acc1),
  ":",
  sub("-.*$", "", acc2)
)]

merged_data[, group := factor(group, levels = c("G0:G0", "G0:G5", "G5:G5"))]

setcolorder(
  merged_data,
  c(
    "sampleComparison", "acc1", "acc2", "group",
    setdiff(names(merged_data), c("sampleComparison", "acc1", "acc2", "group"))
  )
)

fwrite(
  merged_data,
  file.path(outDir, "02.0_pairwise_divergence_BIseq_ONT_merged.tsv"),
  sep = "\t",
  quote = FALSE
)

cat("\nMerged table columns:\n")
print(names(merged_data))

# =========================
# 2) ADD SNP GENETIC DISTANCE
# =========================
mdist_file <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants/cervil-SNPs-ONT_dendro.mdist.gz"
id_file    <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants/cervil-SNPs-ONT_dendro.mdist.id"

merged_data2 <- copy(merged_data)

ids <- fread(id_file, header = FALSE)
samples <- ids$V2
  
mdist <- as.matrix(fread(cmd = paste("zcat", shQuote(mdist_file)), header = FALSE))
rownames(mdist) <- samples
colnames(mdist) <- samples
  
dist_dt <- as.data.table(as.table(mdist))
setnames(dist_dt, c("acc1", "acc2", "genetic_distance"))
  
dist_dt[, acc1 := as.character(acc1)]
dist_dt[, acc2 := as.character(acc2)]
dist_dt <- dist_dt[acc1 != acc2]
  
# same pair key style as merged_data
dist_dt[, pair1 := fifelse(acc1 < acc2, acc1, acc2)]
dist_dt[, pair2 := fifelse(acc1 < acc2, acc2, acc1)]
dist_dt[, sampleComparison := paste(pair1, pair2, sep = ":")]
  
dist_dt <- unique(dist_dt[, .(sampleComparison, genetic_distance)])
  
merged_data2 <- merge(
    merged_data,
    dist_dt,
    by = "sampleComparison",
    all.x = TRUE)
  
fwrite(merged_data2,
    file.path(outDir, "02.0_pairwise_divergence_BIseq_ONT_merged_with_genetic_distance.tsv"),
    sep = "\t",
    quote = FALSE)
  
# =========================
# 3) SUMMARY TABLE
# =========================
value_cols <- grep("^(SNPs-ONT|C-DMR-|CG-DMR-|genetic_distance$)", names(merged_data2), value = TRUE)

summary_dt <- rbindlist(lapply(value_cols, function(v) {
  data.table(
    variable = v,
    n = sum(!is.na(merged_data2[[v]])),
    mean = mean(merged_data2[[v]], na.rm = TRUE),
    median = median(merged_data2[[v]], na.rm = TRUE),
    min = min(merged_data2[[v]], na.rm = TRUE),
    max = max(merged_data2[[v]], na.rm = TRUE)
  )
}))

fwrite(
  summary_dt,
  file.path(outDir, "02.0_pairwise_divergence_BIseq_ONT_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

# =========================
# 4) PLOTS: DMR DIVERGENCE VS SNP DIFFERENCES
# =========================
dmr_cols <- grep("^(C-DMR|CG-DMR)-(BIseq|ONT)$", names(merged_data2), value = TRUE)

plot_counter <- 1
plot_results <- list()

for (y_var in dmr_cols) {
  
  if ("SNPs-ONT" %in% names(merged_data2)) {
    out_file <- sprintf(
      "02.%02d.%s_vs_SNPs-ONT.pdf",
      plot_counter,
      safe_name(y_var)
    )
    
    plot_results[[out_file]] <- create_lm_plot(
      data = merged_data2,
      x_var = "SNPs-ONT",
      y_var = y_var,
      color_by = "group",
      output_file = out_file,
      x_label = "SNP differences / 1,000",
      y_label = paste0(y_var, " differences / 1,000"),
      outPlot = outPlot,
      outDir = outDir,
      color_map = color_map
    )
    
    plot_counter <- plot_counter + 1
  }
  
  if ("genetic_distance" %in% names(merged_data2)) {
    out_file <- sprintf(
      "02.%02d.%s_vs_genetic-distance.pdf",
      plot_counter,
      safe_name(y_var)
    )
    
    plot_results[[out_file]] <- create_lm_plot(
      data = merged_data2,
      x_var = "genetic_distance",
      y_var = y_var,
      color_by = "group",
      output_file = out_file,
      x_label = "SNP genetic distance",
      y_label = paste0(y_var, " differences / 1,000"),
      outPlot = outPlot,
      outDir = outDir,
      color_map = color_map
    )
    
    plot_counter <- plot_counter + 1
  }
}

# =========================
# 5) OPTIONAL: ONT VS BIseq DMR DIVERGENCE
# =========================
for (dmr_class in c("C-DMR", "CG-DMR")) {
  
  x_var <- paste0(dmr_class, "-BIseq")
  y_var <- paste0(dmr_class, "-ONT")
  
  if (all(c(x_var, y_var) %in% names(merged_data2))) {
    
    out_file <- sprintf(
      "02.%02d.%s_ONT_vs_BIseq.pdf",
      plot_counter,
      safe_name(dmr_class)
    )
    
    plot_results[[out_file]] <- create_lm_plot(
      data = merged_data2,
      x_var = x_var,
      y_var = y_var,
      color_by = "group",
      output_file = out_file,
      x_label = paste0(x_var, " differences / 1,000"),
      y_label = paste0(y_var, " differences / 1,000"),
      outPlot = outPlot,
      outDir = outDir,
      color_map = color_map
    )
    
    plot_counter <- plot_counter + 1
  }
}
