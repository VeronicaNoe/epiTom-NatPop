#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================================================
# Paths
# =========================================================
MERGED_DIR <- "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ad_heritability-merged"
OUT_PLOT   <- MERGED_DIR
OUT_TABLE  <- MERGED_DIR

H2_FILE   <- file.path(MERGED_DIR, "00.1_merged-H2_all-markers.tsv")
BETA_FILE <- file.path(MERGED_DIR, "00.2_merged-beta_sig_all-markers.tsv")

C_BED  <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation/C-DMR.bed"
CG_BED <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation/CG-DMR.bed"

# Tables
OUT_H2_FULL_TABLE <- file.path(OUT_TABLE, "review_DMR-GWAS_00.1_merged-H2_genVar_fullUniverse_allH2.tsv")
OUT_STATUS_WIDE   <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_genVar_status_wide.tsv")
OUT_STATUS_COUNTS <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_genVar_status_counts.tsv")
OUT_COMBO_COUNTS  <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_genVar_combo_counts.tsv")
OUT_LINK_COUNTS   <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_genVar_link_status_counts.tsv")
OUT_QTL_DIST_TSV  <- file.path(OUT_TABLE, "review_DMR-GWAS_genVar_mQTL_number_distribution.tsv")
OUT_QTL_TOTAL_TSV <- file.path(OUT_TABLE, "review_DMR-GWAS_total_mQTL_number_genVar.tsv")

# Plots
OUT_COMP_PDF       <- file.path(OUT_PLOT, "review_DMR-GWAS_01_type_composition_genVar.pdf")
OUT_COMP_CLEAN_PDF <- file.path(OUT_PLOT, "review_DMR-GWAS_01.1_link_status_composition_genVar.pdf")
OUT_BETA_ABS_PDF   <- file.path(OUT_PLOT, "review_DMR-GWAS_04_beta_absValues_genVar.pdf")
OUT_BETA_PDF       <- file.path(OUT_PLOT, "review_DMR-GWAS_05_beta_values_genVar.pdf")
OUT_COMBO_BAR_PDF  <- file.path(OUT_PLOT, "review_DMR-GWAS_06_genVar_combo_status_counts.pdf")
OUT_LINK_BAR_PDF   <- file.path(OUT_PLOT, "review_DMR-GWAS_07_genVar_link_status_counts.pdf")
OUT_QTL_PDF        <- file.path(OUT_PLOT, "review_DMR-GWAS_08_QTL_number_distribution_per_DMR_genVar.pdf")
OUT_QTL_TOTAL_PDF  <- file.path(OUT_PLOT, "review_DMR-GWAS_09_total_mQTL_number_genVar.pdf")

marker_levels <- c("SNP", "TIP", "SV")

# =========================================================
# Helpers
# =========================================================
normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^SL2\\.50ch", "", x)
  x <- sub("^ch", "", x, ignore.case = TRUE)
  sprintf("ch%02d", as.integer(x))
}

read_dmr_bed <- function(file, dmr_type) {
  dt <- fread(file, header = FALSE, select = 1:3)
  setnames(dt, c("chr_raw", "start", "end"))
  dt[, chr := normalize_chr(chr_raw)]
  dt[, DMRtype := dmr_type]
  dt[, DMRid := paste0(chr, "_", dmr_type, "_", start)]
  dt[, .(DMRid, chr, start, end, dmr_class = DMRtype)]
}

combine_gen_var_type <- function(SNP, TIP, SV) {
  states <- as.character(c(SNP, TIP, SV))
  states[is.na(states)] <- "nonSig"
  
  if (any(states == "sig")) return("sig")
  
  n_non <- sum(states == "nonSig")
  n_gif <- sum(states == "nonSig-GIF")
  
  if (n_non == 0 && n_gif == 0) return("nonSig")
  if (n_gif > n_non) return("nonSig-GIF")
  return("nonSig")
}

make_status_combo <- function(SNP, TIP, SV) {
  states <- as.character(c(SNP, TIP, SV))
  states[is.na(states)] <- "nonSig"
  paste(sort(unique(states)), collapse = "|")
}

norm_assoc <- function(dt, marker_source) {
  dt <- as.data.table(copy(dt))
  
  nms <- names(dt)
  nms <- gsub("-", "_", nms)
  setnames(dt, old = names(dt), new = nms)
  
  for (nm in c("chrSNP", "start_loci", "end_loci", "meanBeta", "meanSDbeta", "minPvalue",
               "nSNPs", "DMRpos", "chrDMR", "DMR", "startDMR", "lociType", "lociPosition")) {
    if (!nm %in% names(dt)) dt[, (nm) := NA]
  }
  
  dt[, chrSNP       := as.character(chrSNP)]
  dt[, start_loci   := suppressWarnings(as.integer(start_loci))]
  dt[, end_loci     := suppressWarnings(as.integer(end_loci))]
  dt[, meanBeta     := suppressWarnings(as.numeric(meanBeta))]
  dt[, meanSDbeta   := suppressWarnings(as.numeric(meanSDbeta))]
  dt[, minPvalue    := suppressWarnings(as.numeric(minPvalue))]
  dt[, nSNPs        := suppressWarnings(as.integer(nSNPs))]
  dt[, DMRpos       := as.character(DMRpos)]
  dt[, chrDMR       := suppressWarnings(as.integer(chrDMR))]
  dt[, DMR          := as.character(DMR)]
  dt[, startDMR     := suppressWarnings(as.integer(startDMR))]
  dt[, lociType     := as.character(lociType)]
  dt[, lociPosition := as.character(lociPosition)]
  dt[, marker_source := marker_source]
  
  dt[, .(chrSNP, start_loci, end_loci, meanBeta, meanSDbeta, minPvalue, nSNPs,
         DMRpos, chrDMR, DMR, startDMR, lociType, lociPosition, marker_source)]
}

# =========================================================
# 1) Load merged H2 and beta
# =========================================================
h2      <- fread(H2_FILE)
betaVal <- fread(BETA_FILE)
setDT(h2)
setDT(betaVal)

# Fix unnamed V columns if present
vcols <- intersect(paste0("V", 1:6), names(h2))
if (length(vcols) > 0) {
  setnames(
    h2,
    old = vcols,
    new = c("logVariance", "logLik", "delta", "sigmaGen", "sigmaError", "h2")[seq_along(vcols)]
  )
}

for (v in c("sigmaGen", "sigmaError", "h2")) {
  if (v %in% names(h2)) h2[, (v) := suppressWarnings(as.numeric(get(v)))]
}

for (v in c("beta", "betaSD", "pvalue")) {
  if (v %in% names(betaVal)) betaVal[, (v) := suppressWarnings(as.numeric(get(v)))]
}

betaVal[, absValues := abs(beta)]

h2[, marker_set := as.character(marker_set)]
h2[, dmr_class  := as.character(dmr_class)]
h2[, type       := as.character(type)]

betaVal[, marker_set := as.character(marker_set)]
betaVal[, dmr_class  := as.character(dmr_class)]
betaVal[, type       := as.character(type)]

# Build H2 safely from variance components
h2[, H2 := NA_real_]
if (all(c("sigmaGen", "sigmaError") %in% names(h2))) {
  h2[
    !is.na(sigmaGen) & !is.na(sigmaError) & (sigmaGen + sigmaError) > 0,
    H2 := sigmaGen / (sigmaGen + sigmaError)
  ]
}

# =========================================================
# 2) DMR annotation from BEDs
# =========================================================
c_dmr  <- read_dmr_bed(C_BED,  "C-DMR")
cg_dmr <- read_dmr_bed(CG_BED, "CG-DMR")

dmr_annotation <- rbindlist(list(c_dmr, cg_dmr), use.names = TRUE, fill = TRUE)
dmr_annotation <- unique(dmr_annotation, by = "DMRid")

# =========================================================
# 3) Define universe
# =========================================================
# Keep the same rule as before: define the DMR universe from SNP rows.
dmr_universe <- unique(h2[marker_set == "SNP" & !is.na(DMRid), .(DMRid, dmr_class)])

dmr_universe <- merge(
  dmr_universe,
  dmr_annotation[, .(DMRid, chr, start, end)],
  by = "DMRid",
  all.x = TRUE
)

# =========================================================
# 4) Expand SNP-defined universe to SNP/TIP/SV
# =========================================================
dmr_full <- CJ(DMRid = dmr_universe$DMRid, marker_set = marker_levels, unique = TRUE)
dmr_full <- merge(dmr_full, dmr_universe[, .(DMRid, dmr_class)], by = "DMRid", all.x = TRUE)

# =========================================================
# 5) Build full H2 table WITHOUT max collapse
# =========================================================
h2_obs <- copy(h2)[
  DMRid %in% dmr_universe$DMRid &
    marker_set %in% marker_levels,
  .(DMRid, dmr_class, marker_set, marker_type = type, h2, H2)
]

# Recover dmr_class from the universe if missing
h2_obs <- merge(
  h2_obs,
  dmr_universe[, .(DMRid, dmr_class_universe = dmr_class)],
  by = "DMRid",
  all.x = TRUE
)
h2_obs[is.na(dmr_class), dmr_class := dmr_class_universe]
h2_obs[, dmr_class_universe := NULL]

h2_obs[is.na(marker_type), marker_type := "nonSig"]
h2_obs[is.na(H2), H2 := 0]
h2_obs[is.na(h2), h2 := 0]

# Add one zero-filled row only for marker sets that are completely absent
# for a DMR in the input H2 table.
observed_keys <- unique(h2_obs[, .(DMRid, marker_set)])

h2_missing <- dmr_full[
  !observed_keys,
  on = .(DMRid, marker_set)
]
h2_missing[, marker_type := "nonSig"]
h2_missing[, h2 := 0]
h2_missing[, H2 := 0]

h2_full <- rbindlist(
  list(
    h2_obs[, .(DMRid, dmr_class, marker_set, marker_type, h2, H2)],
    h2_missing[, .(DMRid, dmr_class, marker_set, marker_type, h2, H2)]
  ),
  use.names = TRUE,
  fill = TRUE
)

# =========================================================
# 6) Build SNP/TIP/SV status per DMR and collapse status to gen_var
# =========================================================
status_long <- unique(h2_full[, .(DMRid, dmr_class, marker_set, type = marker_type)])
status_long[, type := as.character(type)]
status_long[, marker_set := as.character(marker_set)]
status_long[, dmr_class := as.character(dmr_class)]

status_long[, type_rank := fifelse(
  type == "sig", 3L,
  fifelse(type == "nonSig-GIF", 2L,
          fifelse(type == "nonSig", 1L, 0L))
)]
setorder(status_long, DMRid, marker_set, -type_rank)
status_long <- status_long[, .SD[1], by = .(DMRid, marker_set)]

status_type_wide <- dcast(
  status_long[, .(DMRid, dmr_class, marker_set, type)],
  DMRid + dmr_class ~ marker_set,
  value.var = "type",
  fill = "nonSig"
)

for (x in marker_levels) {
  if (!x %in% names(status_type_wide)) status_type_wide[, (x) := "nonSig"]
}

status_type_wide[, combo_class := mapply(make_status_combo, SNP, TIP, SV)]
status_type_wide[, gen_var_type := mapply(combine_gen_var_type, SNP, TIP, SV)]
status_type_wide[, link_status := fifelse(gen_var_type == "sig", "linked", "nonLinked")]
status_type_wide[, genetic_class := "gen_var"]

# DMR-level table for composition/counts/QTL filtering
h2_gen_status <- copy(status_type_wide)
h2_gen_status[, type := gen_var_type]
h2_gen_status[, dmr_class := factor(dmr_class, levels = c("CG-DMR", "C-DMR"))]
h2_gen_status[, type := factor(type, levels = c("sig", "nonSig", "nonSig-GIF"))]
h2_gen_status[, link_status := factor(link_status, levels = c("linked", "nonLinked"))]
h2_gen_status[, genetic_class := factor(genetic_class, levels = "gen_var")]

# Long H2 table for violin plots: all marker-specific values are retained
h2_gen_all <- merge(
  h2_full,
  status_type_wide[, .(DMRid, dmr_class, SNP, TIP, SV, combo_class, gen_var_type, link_status, genetic_class)],
  by = c("DMRid", "dmr_class"),
  all.x = TRUE
)

h2_gen_all[, type := gen_var_type]
h2_gen_all[, marker_source := marker_set]

h2_gen_all[, dmr_class := factor(dmr_class, levels = c("CG-DMR", "C-DMR"))]
h2_gen_all[, type := factor(type, levels = c("sig", "nonSig"))]
h2_gen_all[, link_status := factor(link_status, levels = c("linked", "nonLinked"))]
h2_gen_all[, genetic_class := factor(genetic_class, levels = "gen_var")]
h2_gen_all[, marker_source := factor(marker_source, levels = marker_levels)]

# Save 
fwrite(h2_gen_all, OUT_H2_FULL_TABLE, sep = "\t")
fwrite(status_type_wide, OUT_STATUS_WIDE, sep = "\t")

status_counts <- h2_gen_status[, .N, by = .(dmr_class, type)][order(dmr_class, -N)]
combo_counts  <- h2_gen_status[, .N, by = .(dmr_class, combo_class)][order(dmr_class, -N)]
link_counts   <- h2_gen_status[, .N, by = .(dmr_class, link_status)][order(dmr_class, -N)]

fwrite(status_counts, OUT_STATUS_COUNTS, sep = "\t")
fwrite(combo_counts,  OUT_COMBO_COUNTS,  sep = "\t")
fwrite(link_counts,   OUT_LINK_COUNTS,   sep = "\t")

print(status_counts)
print(combo_counts)
print(link_counts)

# =========================================================
# 7) Composition plot: gen_var type
# =========================================================
summary_data <- h2_gen_status[!is.na(dmr_class) & !is.na(type), .(count = .N), by = .(dmr_class, genetic_class, type)]
summary_data[, percentage := 100 * count / sum(count), by = .(dmr_class, genetic_class)]

p <- ggplot(summary_data, aes(x = genetic_class, y = percentage, fill = type)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~dmr_class, nrow = 1) +
  scale_fill_manual(values = c(sig = "darkred", nonSig = "#5ba48d", `nonSig-GIF` = "#cccccc")) +
  labs(x = "Genetic variants", y = "Percentage", fill = "Association class") +
  theme_minimal()

ggsave(OUT_COMP_PDF, p, width = 18, height = 12, units = "cm")

# =========================================================
# 7b) Composition plot: linked vs nonLinked
# =========================================================
summary_data_linked <- h2_gen_status[!is.na(dmr_class) & !is.na(link_status), .(count = .N), by = .(dmr_class, genetic_class, link_status)]
summary_data_linked[, percentage := 100 * count / sum(count), by = .(dmr_class, genetic_class)]

p <- ggplot(summary_data_linked, aes(x = genetic_class, y = percentage, fill = link_status)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~dmr_class, nrow = 1) +
  scale_fill_manual(values = c(linked = "darkred", nonLinked = "#5ba48d")) +
  labs(x = "Genetic variants", y = "Percentage", fill = "Link status") +
  theme_minimal()

ggsave(OUT_COMP_CLEAN_PDF, p, width = 18, height = 12, units = "cm")

# =========================================================
# 8) H2 / h2 violin plots using ALL H2/h2 values
# =========================================================
plot_h2 <- copy(h2_gen_all)

p <- ggplot(plot_h2, aes(x = genetic_class, y = H2, fill = dmr_class, color = dmr_class)) +
  geom_violin(alpha = 0.25, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.2, color = "black") +
  scale_fill_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  scale_color_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Genetic variants", y = "H2") +
  facet_grid(dmr_class ~ type) +
  theme_minimal() +
  theme(legend.position = "none")

OUT_H2_PDF         <- file.path(OUT_PLOT, "review_DMR-GWAS_02.1_H2_violin_genVar_allH2.pdf")
ggsave(OUT_H2_PDF, p, width = 24, height = 16, units = "cm")

p <- ggplot(plot_h2, aes(x = genetic_class, y = h2, fill = dmr_class, color = dmr_class)) +
  geom_violin(alpha = 0.25, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.2, color = "black") +
  scale_fill_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  scale_color_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Genetic variants", y = "h2") +
  facet_grid(dmr_class ~ type) +
  theme_minimal() +
  theme(legend.position = "none")

OUT_h2_PDF         <- file.path(OUT_PLOT, "review_DMR-GWAS_02.2_h2_violin_genVar_allH2.pdf")
ggsave(OUT_h2_PDF, p, width = 24, height = 16, units = "cm")

# Linked/nonLinked H2 / h2
plot_h2 <- plot_h2[as.character(type) != "nonSig-GIF"]
p <- ggplot(plot_h2, aes(x = genetic_class, y = H2, fill = dmr_class, color = dmr_class)) +
  geom_violin(alpha = 0.25, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.2, color = "black") +
  scale_fill_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  scale_color_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Genetic variants", y = "H2") +
  facet_grid(dmr_class ~ link_status) +
  theme_minimal() +
  theme(legend.position = "none")
OUT_H2_CLEAN_PDF   <- file.path(OUT_PLOT, "review_DMR-GWAS_03.1_H2_violin_linked_genVar_allH2.pdf")
ggsave(OUT_H2_CLEAN_PDF, p, width = 18, height = 16, units = "cm")

p <- ggplot(plot_h2, aes(x = genetic_class, y = h2, fill = dmr_class, color = dmr_class)) +
  geom_violin(alpha = 0.25, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.2, color = "black") +
  scale_fill_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  scale_color_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Genetic variants", y = "h2") +
  facet_grid(dmr_class ~ link_status) +
  theme_minimal() +
  theme(legend.position = "none")
OUT_h2_CLEAN_PDF   <- file.path(OUT_PLOT, "review_DMR-GWAS_03.2_h2_violin_linked_genVar_allH2.pdf")
ggsave(OUT_h2_CLEAN_PDF, p, width = 18, height = 16, units = "cm")

# Optional diagnostic summary: number of H2 rows per DMR
message("Rows in DMR-level status table: ", nrow(h2_gen_status))
message("Rows in all-H2 violin table: ", nrow(h2_gen_all))
message("Expected approx. 3x DMR rows if one H2 per marker set: ", 3 * nrow(h2_gen_status))

# =========================================================
# 9) Beta plots: all SNP/TIP/SV significant beta rows collapsed as gen_var
# =========================================================
betaVal_use <- betaVal[DMRid %in% dmr_universe$DMRid]
betaVal_use[, genetic_class := "gen_var"]
betaVal_use[, genetic_class := factor(genetic_class, levels = "gen_var")]
betaVal_use[, dmr_class := factor(dmr_class, levels = c("CG-DMR", "C-DMR"))]

p <- ggplot(betaVal_use, aes(x = genetic_class, y = absValues, fill = dmr_class, color = dmr_class)) +
  geom_violin(trim = TRUE, alpha = 0.25) +
  geom_boxplot(width = 0.12, alpha = 0.5, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  scale_color_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Genetic variants", y = "Absolute beta") +
  facet_wrap(~dmr_class, nrow = 1) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(OUT_BETA_ABS_PDF, p, width = 18, height = 12, units = "cm")

p <- ggplot(betaVal_use, aes(x = genetic_class, y = beta, fill = dmr_class, color = dmr_class)) +
  geom_violin(trim = TRUE, alpha = 0.25) +
  geom_boxplot(width = 0.12, alpha = 0.5, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  scale_color_manual(values = c("C-DMR" = "#820a86", "CG-DMR" = "#ffca7b")) +
  labs(x = "Genetic variants", y = "Beta") +
  facet_wrap(~dmr_class, nrow = 1) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(OUT_BETA_PDF, p, width = 18, height = 12, units = "cm")

# =========================================================
# 10) Combo and link-status count plots
# =========================================================
combo_plot <- copy(combo_counts)
combo_plot[, combo_class := factor(combo_class, levels = combo_plot[, unique(combo_class[order(-N)])])]

p <- ggplot(combo_plot, aes(x = N, y = combo_class, fill = dmr_class)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.75, width = 0.7) +
  scale_fill_manual(values = c(`C-DMR` = "#820a86", `CG-DMR` = "#ffca7b")) +
  labs(x = "Number of DMRs", y = "SNP/TIP/SV status combo", fill = "DMR class") +
  theme_minimal()

ggsave(OUT_COMBO_BAR_PDF, p, width = 12, height = 9)

p <- ggplot(link_counts, aes(x = link_status, y = N, fill = dmr_class)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.75, width = 0.7) +
  geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.35, size = 4) +
  scale_fill_manual(values = c(`C-DMR` = "#820a86", `CG-DMR` = "#ffca7b")) +
  labs(x = "Link status", y = "Number of DMRs", fill = "DMR class") +
  theme_minimal()

ggsave(OUT_LINK_BAR_PDF, p, width = 10, height = 7)

# =========================================================
# 11) Association table per marker source, collapsed to gen_var
# =========================================================
files <- list(
  SNP = "/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/01.0_cis-trans_sig-associations.tsv",
  SV  = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results/01.0_cis-trans_sig-associations.tsv",
  TIP = "/mnt/disk2/vibanez/otherAnalysis/11_DMR-GWAS-TIPs/ae_data-analysis/results/01.0_cis-trans_sig-associations.tsv"
)

assoc_list <- lapply(names(files), function(x) norm_assoc(fread(files[[x]]), x))
assoc_dt <- rbindlist(assoc_list, use.names = TRUE, fill = TRUE)

sig_dmrs <- h2_gen_status[type == "sig", as.character(DMRid)]

assoc_dt <- assoc_dt[
  DMRpos %in% dmr_universe$DMRid &
    DMRpos %in% sig_dmrs
]

mqtl_loci <- unique(assoc_dt[, .(
  DMRid = DMRpos,
  dmr_class = DMR,
  marker_source,
  chrSNP,
  start_loci,
  end_loci,
  lociType,
  lociPosition
)])

mqtl_number <- mqtl_loci[, .(mQTL_number = .N), by = .(DMRid, dmr_class)]
mqtl_number[, genetic_class := "gen_var"]

mqtl_number[, class := fifelse(
  mQTL_number == 1, "1",
  fifelse(mQTL_number <= 10, "2-10",
          fifelse(mQTL_number <= 20, "11-20",
                  fifelse(mQTL_number <= 30, "21-30",
                          fifelse(mQTL_number <= 40, "31-40",
                                  fifelse(mQTL_number <= 50, "41-50",
                                          fifelse(mQTL_number <= 100, "51-100", ">100"))))))
)]

class_levels <- c("1", "2-10", "11-20", "21-30", "31-40", "41-50", "51-100", ">100")
mqtl_number[, class := factor(class, levels = class_levels)]
mqtl_number[, dmr_class := factor(dmr_class, levels = c("CG-DMR", "C-DMR"))]
mqtl_number[, genetic_class := factor(genetic_class, levels = "gen_var")]

plot_dt <- mqtl_number[, .(count = .N), by = .(dmr_class, genetic_class, class)]
plot_dt[, proportion := count / sum(count), by = .(dmr_class, genetic_class)]

plot_grid <- CJ(
  dmr_class = levels(h2_gen_status$dmr_class),
  genetic_class = "gen_var",
  class = class_levels,
  unique = TRUE
)
plot_grid[, dmr_class := factor(dmr_class, levels = c("CG-DMR", "C-DMR"))]
plot_grid[, genetic_class := factor(genetic_class, levels = "gen_var")]
plot_grid[, class := factor(class, levels = class_levels)]

plot_dt2 <- merge(plot_grid, plot_dt, by = c("dmr_class", "genetic_class", "class"), all.x = TRUE)
plot_dt2[is.na(count), count := 0L]
plot_dt2[is.na(proportion), proportion := 0]

fwrite(as.data.table(plot_dt2), OUT_QTL_DIST_TSV, sep = "\t")

p <- ggplot(plot_dt2, aes(x = proportion, y = class, fill = dmr_class, color = dmr_class)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.75, width = 0.7) +
  scale_fill_manual(values = c(`C-DMR` = "#820a86", `CG-DMR` = "#ffca7b")) +
  scale_color_manual(values = c(`C-DMR` = "#820a86", `CG-DMR` = "#ffca7b")) +
  labs(x = "Proportion of DMRs", y = "Number of QTLs", fill = "DMR class", color = "DMR class") +
  theme_minimal()

ggsave(OUT_QTL_PDF, p, width = 12, height = 8)

# Total number of mQTLs for gen_var
mqtl_total <- mqtl_number[, .(
  total_mQTL = sum(mQTL_number),
  n_DMR_with_mQTL = uniqueN(DMRid)
), by = .(dmr_class, genetic_class)]

fwrite(as.data.table(mqtl_total), OUT_QTL_TOTAL_TSV, sep = "\t")

p_total_mqtl <- ggplot(mqtl_total, aes(x = genetic_class, y = total_mQTL, fill = dmr_class, color = dmr_class)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.75, width = 0.7) +
  geom_text(aes(label = format(total_mQTL, big.mark = ",")), vjust = -0.35, size = 4) +
  scale_fill_manual(values = c(`C-DMR` = "#820a86", `CG-DMR` = "#ffca7b")) +
  scale_color_manual(values = c(`C-DMR` = "#820a86", `CG-DMR` = "#ffca7b")) +
  labs(x = "Genetic variants", y = "Total number of mQTLs", fill = "DMR class", color = "DMR class") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11)
  )

ggsave(OUT_QTL_TOTAL_PDF, p_total_mqtl, width = 10, height = 7)
