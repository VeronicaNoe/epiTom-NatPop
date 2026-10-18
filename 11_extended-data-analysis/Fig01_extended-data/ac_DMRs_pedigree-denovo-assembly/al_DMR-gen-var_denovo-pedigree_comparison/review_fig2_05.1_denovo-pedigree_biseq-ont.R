#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
genvar <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants/results/02.05_variable_variants.all.genotypes.tsv"

ONT_data <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ai_DMR-classification/ba_get-DMR/"
BIseq_data <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.0_get-DMR/ab_predigree/cervil-denovo/"

# Use the directory, not the single C-DMR file
Biseq_sl25_denovo_dir <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation/cervil-sl25/"

outDir <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/al_DMR-gen-var_denovo-pedigree_comparison/"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
standardize_dmr <- function(x) {
  x <- gsub("\\.bed$", "", x)
  x <- gsub("\\.gz$", "", x)
  
  if (x %in% c("C", "C-DMR")) return("C-DMR")
  if (x %in% c("CG", "CG-DMR")) return("CG-DMR")
  
  return(NA_character_)
}

parse_file_info <- function(file_name) {
  base <- basename(file_name)
  base <- gsub("\\.bed$", "", base)
  base <- gsub("\\.gz$", "", base)
  
  parts <- strsplit(base, "_", fixed = TRUE)[[1]]
  
  if (length(parts) < 3) return(NULL)
  
  p1 <- strsplit(parts[1], "-", fixed = TRUE)[[1]]
  p2 <- strsplit(parts[2], "-", fixed = TRUE)[[1]]
  
  if (length(p1) < 2 || length(p2) < 2) return(NULL)
  
  g1 <- p1[1]
  s1 <- p1[2]
  g2 <- p2[1]
  s2 <- p2[2]
  
  dmr_class <- standardize_dmr(parts[3])
  
  if (is.na(dmr_class)) return(NULL)
  
  comparison <- ifelse(g1 == g2, paste0("intra-", g1), "inter-Gen")
  
  list(
    g1 = g1,
    s1 = s1,
    g2 = g2,
    s2 = s2,
    DMR = dmr_class,
    comparison = comparison,
    samplePair = paste0(g1, ":", s1, "-", g2, ":", s2),
    sample = paste0(s1, "-", s2)
  )
}

read_lifted_dmr <- function(dmr_class, lifted_dir) {
  f <- file.path(lifted_dir, paste0(dmr_class, ".SL25_to_polished.bed"))
  
  if (!file.exists(f)) {
    stop("Missing lifted DMR file: ", f)
  }
  
  dt <- fread(
    f,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    data.table = TRUE
  )
  
  if (ncol(dt) < 4) {
    stop("Lifted DMR file has fewer than 4 columns: ", f)
  }
  
  dt <- dt[, .(
    DMR = dmr_class,
    chr = as.character(V1),
    start = as.integer(V2),
    end = as.integer(V3),
    dmr_id = as.character(V4)
  )]
  
  dt <- unique(dt)
  dt
}


read_pairwise_dmr_dir <- function(inDir, dataset_name, lifted_dt) {
  input <- list.files(path = inDir,    pattern = "\\.bed$|\\.bed\\.gz$",    full.names = TRUE  )
  df_list <- list()
  for (i in input) {
    info <- parse_file_info(i)
    tmp <- tryCatch({
      fread(
        i,
        sep = "\t",
        header = FALSE,
        data.table = TRUE,
        fill = TRUE,
        na.strings = c(""),
        nThread = 20)})
    binTmp <- tmp[, .(
      chr = as.character(V1),
      start = as.integer(V2),
      end = as.integer(V3),
      pvalue = as.numeric(V11),
      qvalue = as.numeric(V12),
      difMeth = as.numeric(V13),
      sample1 = fifelse(as.numeric(V5) > 0, as.numeric(V6) / as.numeric(V5), NA_real_),
      sample2 = fifelse(as.numeric(V8) > 0, as.numeric(V9) / as.numeric(V8), NA_real_))]
    
    binTmp[, `:=`(
      dataset = dataset_name,
      file = basename(i),
      samplePair = info$samplePair,
      DMR = info$DMR,
      comparison = info$comparison,
      sample = info$sample)]
    
    df_list[[length(df_list) + 1]] <- binTmp
  }
  df <- rbindlist(df_list, fill = TRUE)
  # Significant DMRs only
  #df <- df[qvalue <= 0.005]
  df[, pos := paste(chr, start, end, sep = "-")]
  setorder(df, dataset, DMR, chr, start, end)
  
  df
}

# -----------------------------
# Read lifted C-DMR and CG-DMR reference intervals
# -----------------------------
lifted_dt <- rbindlist(list(read_lifted_dmr("C-DMR", Biseq_sl25_denovo_dir), 
                            read_lifted_dmr("CG-DMR", Biseq_sl25_denovo_dir)), fill = TRUE)

# -----------------------------
# Read ONT and BIseq pairwise DMRs
# -----------------------------
df_ont <- read_pairwise_dmr_dir(inDir = ONT_data, dataset_name = "ONT",  lifted_dt = lifted_dt)
# Keep only DMRs overlapping the lifted SL25 -> polished BIseq DMR set

df_ont[, row_id := .I]
setkey(lifted_dt, DMR, chr, start, end)
setkey(df_ont, DMR, chr, start, end)
hits <- foverlaps(df_ont,lifted_dt,
                  by.x = c("DMR", "chr", "start", "end"),
                  by.y = c("DMR", "chr", "start", "end"),
                  type = "any",nomatch = 0L)
hit_ids <- hits[, .(dmr_id = paste(unique(dmr_id), collapse = ";")), by = row_id]
df_ont_clean <- merge(df_ont, hit_ids, by = "row_id")
df_ont_clean[, row_id := NULL]

df_biseq <- read_pairwise_dmr_dir(inDir = BIseq_data,  dataset_name = "BIseq",  lifted_dt = lifted_dt)
df <- rbindlist(list(df_ont_clean, df_biseq), fill = TRUE)

fwrite(out,file = file.path(outDir, "01.0_pairwise_DMRs_ONT_BIseq.tsv"),sep = "\t",  quote = FALSE)

# -----------------------------
# Pairwise DMR count summary
# -----------------------------
pairWiseStatSumm <- df %>%
  group_by(dataset, DMR, comparison, samplePair) %>%
  summarise(nDMR = n(), .groups = "drop") %>%
  mutate(generations_numeric = ifelse(comparison == "inter-Gen", 6, 2))


# Add generation 0 baseline for intra-generation comparisons
pairWiseStatSumm <- pairWiseStatSumm %>%
  bind_rows(pairWiseStatSumm %>%
      filter(generations_numeric == 2) %>%
      mutate(generations_numeric = 0, nDMR = 0)  )

# Regression-like divergence rates
safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

lm_results <- pairWiseStatSumm %>%
  group_by(dataset, DMR) %>%
  summarise(
    mean_g2 = safe_mean(nDMR[generations_numeric == 2]),
    mean_g6 = safe_mean(nDMR[generations_numeric == 6]),
    beta_early = mean_g2 / 2,
    beta_late = (mean_g6 - mean_g2) / 4,
    max_y = max(nDMR, na.rm = TRUE),
    y_position_early = pmax(max_y * 0.85, 1),
    y_position_late = pmax(max_y * 0.70, 1),
    .groups = "drop"  )

pairWiseStatSumm <- left_join(pairWiseStatSumm,  lm_results,  by = c("dataset", "DMR"))
fwrite(pairWiseStatSumm,file = file.path(outDir, "01.1_number-DMR-pairwise_ONT_BIseq.filtered_by_SL25_lifted_DMRs.tsv"),
       sep = "\t",  quote = FALSE)

# -----------------------------
# Plot
# -----------------------------
dmr_colors <- c( "C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")

p <- ggplot( pairWiseStatSumm,aes(x = generations_numeric, y = nDMR, color = DMR)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(data = subset(pairWiseStatSumm, generations_numeric <= 2),
    method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", linewidth = 0.7) +
  geom_smooth(data = subset(pairWiseStatSumm, generations_numeric >= 2),
    method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", linewidth = 0.7) +
  facet_grid(DMR ~ dataset, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 6, by = 2), limits = c(0, 6)) +
  scale_color_manual(values = dmr_colors) +
  labs(title = "# of DMRs by pairwise comparison", x = "Generations",    y = "# of DMRs"  ) +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  geom_text(data = lm_results, aes(x = 1, y = y_position_early, label = paste0("DMR div/gen = ", round(beta_early, 3))),
    inherit.aes = FALSE, size = 4, fontface = "bold", color = "black") +
  geom_text(data = lm_results, aes(x = 4, y = y_position_late, label = paste0("DMR div/gen = ", round(beta_late, 3)) ),
            inherit.aes = FALSE, size = 4, fontface = "bold",  color = "black"  )

ggsave(filename = file.path(outDir, "01.1_DMR-divergence-per-generation_ONT_BIseq.pdf"),plot = p,  width = 10,  height = 7)
##
# -----------------------------
# Genetic variant pairwise divergence inside lifted DMRs
# -----------------------------
genvar <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants/results/02.05_variable_variants.all.genotypes.tsv"
samples <- c("G0-C3", "G0-C6", "G0-C9", "G5-C1", "G5-C2", "G5-C3")
gt_cols <- paste0(samples, "_normGT")

gv <- fread(genvar, sep = "\t", header = TRUE, data.table = TRUE)

# Standardize coordinates
gv <- gv[, .(
  chr = as.character(CHROM),
  start = as.integer(start),
  end = as.integer(end),
  variant_id,
  variant_class,
  variant_type,
  G0_C3 = `G0-C3_normGT`,
  G0_C6 = `G0-C6_normGT`,
  G0_C9 = `G0-C9_normGT`,
  G5_C1 = `G5-C1_normGT`,
  G5_C2 = `G5-C2_normGT`,
  G5_C3 = `G5-C3_normGT`
)]

setnames(gv,  old = c("G0_C3", "G0_C6", "G0_C9", "G5_C1", "G5_C2", "G5_C3"),  new = samples)

# -----------------------------
# Create pairwise genotype differences
# -----------------------------
is_missing_gt <- function(x) {
  is.na(x) | x %in% c("", ".", "./.", ".|.")
}

pairs <- combn(samples, 2, simplify = FALSE)
genvar_pairwise_list <- list()
for (p in pairs) {
  
  s1 <- p[1]
  s2 <- p[2]
  
  g1 <- strsplit(s1, "-", fixed = TRUE)[[1]][1]
  g2 <- strsplit(s2, "-", fixed = TRUE)[[1]][1]
  
  c1 <- s1
  c2 <- s2
  
  comparison <- ifelse(g1 == g2, paste0("intra-", g1), "inter-Gen")
  samplePair <- paste0(g1, ":", c1, "-", g2, ":", c2)
  sample_name <- paste0(c1, "-", c2)
  
  tmp <- copy(gv)
  
  tmp[, gt1 := get(s1)]
  tmp[, gt2 := get(s2)]
  
  tmp <- tmp[
    !is_missing_gt(gt1) &
      !is_missing_gt(gt2) &
      gt1 != gt2
  ]
  
  if (nrow(tmp) == 0) next
  
  tmp[, `:=`(
    dataset = "GenVar",
    samplePair = samplePair,
    comparison = comparison,
    sample = sample_name,
    sample1 = s1,
    sample2 = s2
  )]
  
  genvar_pairwise_list[[length(genvar_pairwise_list) + 1]] <- tmp
}

genvar_pairwise <- rbindlist(genvar_pairwise_list, fill = TRUE)

# -----------------------------
# Summary for plot
# -----------------------------
pairWiseStatSumm_genvar <- genvar_pairwise %>%
  group_by(dataset,comparison, samplePair) %>%
  summarise(nVar = uniqueN(variant_id), .groups = "drop") %>%
  mutate(generations_numeric = ifelse(comparison == "inter-Gen", 6, 2))

# Add generation 0 baseline
pairWiseStatSumm_genvar <- pairWiseStatSumm_genvar %>%
  bind_rows(pairWiseStatSumm_genvar %>%
              filter(generations_numeric == 2) %>%
              mutate(generations_numeric = 0, nVar = 0)  )

# -----------------------------
# Divergence rates
# -----------------------------
safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

lm_results_genvar <- pairWiseStatSumm_genvar %>%
  group_by(dataset) %>%
  summarise(
    mean_g2 = safe_mean(nVar[generations_numeric == 2]),
    mean_g6 = safe_mean(nVar[generations_numeric == 6]),
    beta_early = mean_g2 / 2,
    beta_late = (mean_g6 - mean_g2) / 4,
    beta_mean = mean_g6 / 6,
    beta_lm = coef(lm(nVar ~ generations_numeric))[2],
    max_y = max(nVar, na.rm = TRUE),
    y_position = pmax(max_y * 0.9, 1),
    y_position_early = pmax(max_y * 0.85, 1),
    y_position_late = pmax(max_y * 0.70, 1),
    .groups = "drop"
  )

pairWiseStatSumm_genvar <- left_join(pairWiseStatSumm_genvar,  lm_results_genvar, by = c("dataset"))
fwrite(pairWiseStatSumm_genvar,  file = file.path(outDir, "01.1_number-genvar-pairwise_inside_lifted_DMRs.tsv"),  sep = "\t",  quote = FALSE)

# -----------------------------
# Plot
# -----------------------------
p_genvar <- ggplot( pairWiseStatSumm_genvar, aes(x = generations_numeric, y = nVar, color = dataset)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth( method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", linewidth = 0.7) +
  scale_x_continuous(breaks = seq(0, 6, by = 2), limits = c(0, 6)) +
  labs(title = "# of genetic variant differences by pairwise comparison",subtitle = "SNPs variants across pedigree samples",
    x = "Generations", y = "# of divergent genetic variants") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  geom_text(data = lm_results_genvar, aes(x = 3, y = y_position, label = paste0("GenVar div/gen = ", round(beta_lm, 3))),
    inherit.aes = FALSE, size = 4, fontface = "bold", color = "black")

ggsave(filename = file.path(outDir, "01.1_genvar-divergence-per-generation.pdf"), plot = p_genvar,  width = 8,  height = 6)

# 
p_genvar <- ggplot(pairWiseStatSumm_genvar,  aes(x = generations_numeric, y = nVar, color = dataset)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(data = subset(pairWiseStatSumm_genvar, generations_numeric <= 2),
              method = "lm",
              formula = y ~ x,
              se = FALSE,
              linetype = "dashed",
              linewidth = 0.7  ) +
  geom_smooth( data = subset(pairWiseStatSumm_genvar, generations_numeric >= 2),
               method = "lm",
               formula = y ~ x,
               se = FALSE,
               linetype = "dashed",
               linewidth = 0.7  ) +
  scale_x_continuous( breaks = seq(0, 6, by = 2),    limits = c(0, 6)) +
  labs(title = "# of genetic variant differences by pairwise comparison",
       subtitle = "SNPs",
       x = "Generations",    y = "# of divergent genetic variants"  ) +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold")  ) +
  geom_text(data = lm_results_genvar,
            aes(x = 1,      y = y_position_early,      label = paste0("var div/gen = ", round(beta_early, 3))    ),
            inherit.aes = FALSE,    size = 4,    fontface = "bold",    color = "black"  ) +
  geom_text(data = lm_results_genvar, aes(x = 4, y = y_position_late, label = paste0("var div/gen = ", round(beta_late, 3))    ),
            inherit.aes = FALSE,
            size = 4,
            fontface = "bold",
            color = "black")

ggsave(filename = file.path(outDir, "01.1_genvar-divergence_by_points.pdf"),  plot = p_genvar,  width = 8,  height = 6)
