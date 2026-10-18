#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ComplexUpset)
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

OUT_H2_FULL_TABLE   <- file.path(OUT_TABLE, "review_DMR-GWAS_00.1_merged-H2_all-markers_fullUniverse_SNPuniverse.tsv")
OUT_UPSET_TABLE     <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_universe_with_link_status_allMarkers_SNPuniverse.tsv")
OUT_UPSET_COUNTS    <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_link_status_counts_allMarkers_SNPuniverse.tsv")
OUT_STATUS_WIDE     <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_marker_type_status_wide_SNPuniverse.tsv")
OUT_STATUS_COUNTS   <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_marker_type_status_counts_SNPuniverse.tsv")
OUT_STATUS_COMBO    <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_marker_type_status_combo_counts_SNPuniverse.tsv")
OUT_STATUS_UPSET    <- file.path(OUT_TABLE, "review_DMR-GWAS_DMR_marker_sig_binary_wide_SNPuniverse.tsv")
OUT_QTL_TOTAL_TSV <- file.path(OUT_TABLE,"review_DMR-GWAS_total_mQTL_number_per_marker_allMarkers_SNPuniverse.tsv")

OUT_COMP_PDF        <- file.path(OUT_PLOT,  "review_DMR-GWAS_01_type_composition_all-markers_SNPuniverse.pdf")
OUT_H2_PDF          <- file.path(OUT_PLOT,  "review_DMR-GWAS_02_H2_violin_all-markers_SNPuniverse.pdf")
OUT_h2_PDF          <- file.path(OUT_PLOT,  "review_DMR-GWAS_03_h2_violin_all-markers_SNPuniverse.pdf")
OUT_BETA_ABS_PDF    <- file.path(OUT_PLOT,  "review_DMR-GWAS_04_beta_absValues_all-markers_SNPuniverse.pdf")
OUT_BETA_PDF        <- file.path(OUT_PLOT,  "review_DMR-GWAS_05_beta_values_all-markers_SNPuniverse.pdf")
OUT_UPSET_PDF       <- file.path(OUT_PLOT,  "review_DMR-GWAS_06_UpSet_DMR_retained_SNP_SV_TIP_allMarkers_SNPuniverse.pdf")
OUT_UPSET_CLEAN_PDF       <- file.path(OUT_PLOT,  "review_DMR-GWAS_07_UpSet_DMR_retained_SNP_SV_TIP_allMarkers_SNPuniverse.pdf")
OUT_QTL_PDF         <- file.path(OUT_PLOT,  "review_DMR-GWAS_08_QTL_number_distribution_per_DMR_allMarkers_SNPuniverse.pdf")
OUT_QTL_TOTAL_PDF <- file.path( OUT_PLOT,"review_DMR-GWAS_09_total_mQTL_number_per_marker_allMarkers_SNPuniverse.pdf")

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
  dt[, .(DMRid, chr, start, end, DMRtype)]
}

make_combo <- function(SNP, SV, TIP, Retained) {
  if (Retained == 1) return("Retained_only")
  labs <- c("SNP", "SV", "TIP")[c(SNP, SV, TIP) == 1]
  paste(labs, collapse = "+")
}

combine_dmr_type <- function(snp, tip, sv) {
  states <- c(snp, tip, sv)
  
  if (any(states == "sig")) return("sig")
  
  n_non    <- sum(states == "nonSig")
  n_gif    <- sum(states == "nonSig-GIF")
  n_absent <- sum(states == "absent")
  
  if (n_non == 0 && n_gif == 0 && n_absent == 3) return("absent")
  if (n_non == 0 && n_gif == 0) return("retained")
  
  if (n_gif > n_non) return("nonSig-GIF")
  if (n_non > n_gif) return("nonSig")
  
  return("nonSig-GIF")
}

make_upset <- function(dt, retained_col, title_text) {
  upset(dt,
    intersect = c("Retained", "All_nonSig_GIF", "SV", "TIP", "SNP"),
    name = "DMRs",
    min_size = 0,
    sort_intersections_by = "cardinality",
    sort_sets = FALSE,
    base_annotations = list("Intersection size" = intersection_size(text = list(size = 3))),
    set_sizes = upset_set_size() + xlab("Set size"),
    queries = list(
      upset_query(intersect = "Retained", fill = retained_col, color = retained_col),
      upset_query(intersect = "SNP", fill = "#999999", color = "#999999"),
      upset_query(intersect = "TIP", fill = "black", color = "black"),
      upset_query(intersect = "SV",  fill = "#977777", color = "#977777"),
      upset_query(intersect = c("SV", "TIP"), fill = "#4D4D4D", color = "#4D4D4D"),
      upset_query(intersect = c("SNP", "SV"), fill = "#C9A3A3", color = "#C9A3A3"),
      upset_query(intersect = c("SNP", "TIP"), fill ="#4D4D4D", color = "#4D4D4D"),
      upset_query(intersect = c("SNP", "SV", "TIP"), fill = "#994301", color = "#994301"),
      upset_query(intersect = "All_nonSig_GIF", fill = "white", color = "black"))) +
    ggtitle(title_text) +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(size = 11) )
}

make_upset_clean <- function(dt, retained_col, title_text) {
  upset(
    dt,
    intersect = c("Retained", "SV", "TIP", "SNP"),
    name = "DMRs",
    min_size = 0,
    sort_intersections_by = "cardinality",
    sort_sets = FALSE,
    base_annotations = list(
      "Intersection size" = intersection_size(text = list(size = 3))
    ),
    set_sizes = upset_set_size() + xlab("Set size"),
    queries = list(
      upset_query(intersect = "Retained", fill = retained_col, color = retained_col),
      upset_query(intersect = "SNP", fill = "#999999", color = "#999999"),
      upset_query(intersect = "TIP", fill = "black", color = "black"),
      upset_query(intersect = "SV",  fill = "#977777", color = "#977777"),
      upset_query(intersect = c("SV", "TIP"), fill = "#4D4D4D", color = "#4D4D4D"),
      upset_query(intersect = c("SNP", "SV"), fill = "#C9A3A3", color = "#C9A3A3"),
      upset_query(intersect = c("SNP", "TIP"), fill = "#4D4D4D", color = "#4D4D4D"),
      upset_query(intersect = c("SNP", "SV", "TIP"), fill = "#994301", color = "#994301")
    )
  ) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 11)
    )
}

# =========================================================
# 1) Load merged H2 and beta
# =========================================================
h2 <- fread(H2_FILE)
betaVal <- fread(BETA_FILE)

setnames(h2, old = intersect(c("V1","V2","V3","V4","V5","V6"), names(h2)),new = c("logVariance","logLik","delta","sigmaGen","sigmaError","h2")[seq_len(length(intersect(c("V1","V2","V3","V4","V5","V6"), names(h2))))])

setDT(h2)
setDT(betaVal)

# numeric cleanup
for (v in c("sigmaGen", "sigmaError", "h2")) {
  if (v %in% names(h2)) h2[, (v) := suppressWarnings(as.numeric(get(v)))]
}

for (v in c("beta", "betaSD", "pvalue")) {
  if (v %in% names(betaVal)) betaVal[, (v) := suppressWarnings(as.numeric(get(v)))]
}

betaVal[, absValues := abs(beta)]

# use characters first to avoid factor issues during joins
h2[, marker_set := as.character(marker_set)]
h2[, dmr_class  := as.character(dmr_class)]
h2[, type       := as.character(type)]

betaVal[, marker_set := as.character(marker_set)]
betaVal[, dmr_class  := as.character(dmr_class)]
betaVal[, type       := as.character(type)]

# build H2 safely
h2[, H2 := NA_real_]
if (all(c("sigmaGen", "sigmaError") %in% names(h2))) {
  h2[
    !is.na(sigmaGen) & !is.na(sigmaError) & (sigmaGen + sigmaError) > 0,
    H2 := sigmaGen / (sigmaGen + sigmaError)
  ]
}

# =========================================================
# 2) DMR annotation from BEDs
#    only for coordinates/annotation, NOT for universe definition
# =========================================================
c_dmr  <- read_dmr_bed(C_BED,  "C-DMR")
cg_dmr <- read_dmr_bed(CG_BED, "CG-DMR")

dmr_annotation <- rbindlist(list(c_dmr, cg_dmr), use.names = TRUE, fill = TRUE)
dmr_annotation <- unique(dmr_annotation, by = "DMRid")
dmr_annotation <- dmr_annotation[, .(DMRid, chr, start, end, dmr_class = DMRtype)]

# =========================================================
# 3) Define universe from SNP DMR IDs only
# =========================================================
dmr_universe <- unique(h2[marker_set == "SNP" & !is.na(DMRid), .(DMRid, dmr_class)])

# optional coordinates from BED annotation
dmr_universe <- merge(dmr_universe,
                      dmr_annotation[, .(DMRid, chr, start, end)],
                      by = "DMRid",  all.x = TRUE)

# =========================================================
# 4) Expand SNP-defined universe to all marker sets
#    one row per SNP-universe DMRid x marker_set
# =========================================================
dmr_full <- CJ(DMRid = dmr_universe$DMRid, marker_set = marker_levels, unique = TRUE)

dmr_full <- merge(dmr_full,  dmr_universe,  by = "DMRid",  all.x = TRUE)

# =========================================================
# 5) Expand H2 to full SNP universe
# =========================================================
h2_small <- copy(h2)[, .(DMRid, dmr_class, marker_set, type, h2, H2)]

# keep strongest class if duplicated:
# sig > nonSig-GIF > nonSig > absent
h2_small[, type_rank := fifelse(
  type == "sig", 3L,
  fifelse(type == "nonSig-GIF", 2L,
          fifelse(type == "nonSig", 1L,
                  fifelse(type == "absent", 0L, -1L)))
)]

setorder(h2_small, DMRid, marker_set, -type_rank,-h2, -H2)
h2_small <- h2_small[, .SD[1], by = .(DMRid, marker_set)]

h2_full <- merge(dmr_full, h2_small[, .(DMRid, marker_set, type,h2, H2)],by = c("DMRid", "marker_set"),  all.x = TRUE)

# within the SNP universe:
h2_full[is.na(type), type := "nonSig"]
h2_full[is.na(H2),   H2   := 0]
h2_full[is.na(h2),   h2   := 0]

h2_full[, marker_set := factor(marker_set, levels = c("SNP", "TIP", "SV"))]
h2_full[, dmr_class  := factor(dmr_class,  levels = c("CG-DMR", "C-DMR"))]
h2_full[, type       := factor(type,       levels = c("sig", "nonSig", "nonSig-GIF"))]

fwrite(h2_full, OUT_H2_FULL_TABLE, sep = "\t")

# =========================================================
# 6) Composition plot
#    based on SNP-defined universe
# =========================================================
summary_data <- h2_full[!is.na(dmr_class) & !is.na(type), .(count = .N),
  by = .(dmr_class, marker_set, type)][, percentage := 100 * count / sum(count),
  by = .(dmr_class, marker_set)]

p <- ggplot(summary_data, aes(x = marker_set, y = percentage, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~dmr_class, nrow = 1) +
  scale_fill_manual(values = c(sig = "darkred", nonSig = "#5ba48d",`nonSig-GIF` = "#cccccc")) +
  labs(x = "Marker set", y = "Percentage", fill = "Association class") +
  theme_minimal()

ggsave(OUT_COMP_PDF, p, width = 18, height = 12, units = "cm")

# =========================================================
# 7) H2 plots
# =========================================================
plot_h2 <- copy(h2_full)
p <- ggplot(plot_h2, aes(x = marker_set, y = H2, fill = marker_set, color = marker_set)) +
  geom_violin(alpha = 0.25, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.2, color = "black") +
  scale_fill_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  scale_color_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Marker set", y = "H2") +
  facet_grid(dmr_class ~ type) +
  theme_minimal()

ggsave(OUT_H2_PDF, p, width = 28, height = 16, units = "cm")

p <- ggplot(plot_h2, aes(x = marker_set, y = h2, fill = marker_set, color = marker_set)) +
  geom_violin(alpha = 0.25, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.2, color = "black") +
  scale_fill_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  scale_color_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Marker set", y = "H2") +
  facet_grid(dmr_class ~ type) +
  theme_minimal()

ggsave(OUT_h2_PDF, p, width = 28, height = 16, units = "cm")
# =========================================================
# 8) Beta plots
#    restrict beta to SNP-universe DMRs for fair comparison
# =========================================================
betaVal_use <- betaVal[DMRid %in% dmr_universe$DMRid]

betaVal_use[, marker_set := factor(marker_set, levels = c("SNP", "TIP", "SV"))]
betaVal_use[, dmr_class  := factor(dmr_class,  levels = c("CG-DMR", "C-DMR"))]

p <- ggplot(betaVal_use, aes(x = marker_set, y = absValues, fill = marker_set, color = marker_set)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  geom_boxplot(width = 0.12, alpha = 0.5, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  scale_color_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  labs(x = "Marker set", y = "Absolute beta") +
  facet_wrap(~dmr_class, nrow = 1) +
  theme_minimal()

ggsave(OUT_BETA_ABS_PDF, p, width = 18, height = 12, units = "cm")

p <- ggplot(betaVal_use, aes(x = marker_set, y = beta, fill = marker_set, color = marker_set)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  geom_boxplot(width = 0.12, alpha = 0.5, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  scale_color_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  labs(x = "Marker set", y = "Beta") +
  facet_wrap(~dmr_class, nrow = 1) +
  theme_minimal()

ggsave(OUT_BETA_PDF, p, width = 18, height = 12, units = "cm")

# =========================================================
# 9) Build all-marker status per DMR
#    from SNP-defined universe
# =========================================================
status_long <- unique(h2_full[, .(DMRid, dmr_class, marker_set, type)])

status_marker_counts <- status_long[, .N, by = .(marker_set, dmr_class, type)][
  order(marker_set, dmr_class, type)]
print(status_marker_counts)

status_long[, type := as.character(type)]
status_long[, marker_set := as.character(marker_set)]
status_long[, dmr_class := as.character(dmr_class)]

status_long[, type_rank := fifelse(type == "sig", 3L,
                                   fifelse(type == "nonSig-GIF", 2L,
                                           fifelse(type == "nonSig", 1L, 0L)))]

setorder(status_long, DMRid, marker_set, -type_rank)
status_long <- status_long[, .SD[1], by = .(DMRid, marker_set)]

status_type_wide <- dcast(status_long[, .(DMRid, dmr_class, marker_set, type)],
  DMRid + dmr_class ~ marker_set, value.var = "type",  fill = "absent") #check absent -> not necessary

for (x in c("SNP", "TIP", "SV")) {
  if (!x %in% names(status_type_wide)) status_type_wide[, (x) := "absent"]
}

status_type_wide[, combined_type := mapply(combine_dmr_type, SNP, TIP, SV)]
#

combo_counts <- status_type_wide %>%
  rowwise() %>%
  mutate(combo = paste(sort(unique(c(SNP, SV, TIP))), collapse = "|")) %>%
  ungroup() %>%
  count(dmr_class, combo, sort = TRUE)

combo_counts
#
fwrite(status_type_wide, OUT_STATUS_WIDE, sep = "\t")

status_combined_counts <- status_type_wide[, .N, by = .(dmr_class, combined_type)][
  order(dmr_class, -N)
]
fwrite(status_combined_counts, OUT_STATUS_COUNTS, sep = "\t")
print(status_combined_counts)

status_combo_counts <- status_type_wide[, .N, by = .(dmr_class, SNP, TIP, SV, combined_type)][order(dmr_class, -N)]
fwrite(status_combo_counts, OUT_STATUS_COMBO, sep = "\t")

# =========================================================
# 10-12) Final UpSet input
#       merge SNP universe first, then classify
#       linked = sig only
# =========================================================
status_base <- copy(status_type_wide)[, .(DMRid,dmr_class, SNP_status = as.character(SNP), TIP_status = as.character(TIP),SV_status  = as.character(SV))]
# Merge with DMR universe BEFORE classification
dmr_all <- merge(dmr_universe,status_base,by = c("DMRid", "dmr_class"),  all.x = TRUE)

# Count status types across marker sets
dmr_all[, n_sig := (SNP_status == "sig") + (TIP_status == "sig") + (SV_status  == "sig")]

dmr_all[, n_gif := (SNP_status == "nonSig-GIF") + (TIP_status == "nonSig-GIF") + (SV_status  == "nonSig-GIF")]

dmr_all[, n_non := (SNP_status == "nonSig") + (TIP_status == "nonSig") + (SV_status  == "nonSig")]

# Binary columns for UpSet: significant links only
dmr_all[, SNP := as.integer(SNP_status == "sig")]
dmr_all[, TIP := as.integer(TIP_status == "sig")]
dmr_all[, SV  := as.integer(SV_status  == "sig")]

# Special non-linked classes
dmr_all[, All_nonSig_GIF := as.integer(n_sig == 0 & n_gif == 3)]
dmr_all[, Retained := as.integer(n_sig == 0 & n_gif < 3)]

# Final mutually interpretable class label
dmr_all[, class_combo := fifelse(n_sig > 0,mapply(make_combo, SNP, SV, TIP, 0),
  fifelse(All_nonSig_GIF == 1,"All_nonSig-GIF","Retained" ))]
# Save status table used for UpSet
fwrite(dmr_all, OUT_UPSET_TABLE, sep = "\t")

# Counts
counts_all <- dmr_all[, .N, by = .(dmr_class, class_combo)][
  order(dmr_class, -N)]

fwrite(counts_all, OUT_UPSET_COUNTS, sep = "\t")
print(counts_all)
# check
print(dmr_all[, .N, by = .(dmr_class, SNP_status,    TIP_status,    SV_status,    class_combo)][order(dmr_class, class_combo, -N)])
# =========================================================
# 13) UpSet PDF
# =========================================================
pdf(OUT_UPSET_PDF, width = 12, height = 7, onefile = TRUE)
print(make_upset( dmr_all[dmr_class == "C-DMR"], retained_col = "#820a86",
                  title_text = "C-DMRs: Retained vs linked to SNP / SV / TIP"))
print(make_upset(dmr_all[dmr_class == "CG-DMR"], retained_col = "#ffca7b",
                 title_text = "CG-DMRs: Retained vs linked to SNP / SV / TIP"))
dev.off()

# clean version
dmr_all_clean <-dmr_all %>%
  filter(class_combo!= "All_nonSig-GIF")
pdf(OUT_UPSET_CLEAN_PDF, width = 12, height = 7, onefile = TRUE)
print(make_upset_clean( dmr_all_clean[dmr_class == "C-DMR"], retained_col = "#820a86",
                  title_text = "C-DMRs: Retained vs linked to SNP / SV / TIP"))
print(make_upset_clean(dmr_all_clean[dmr_class == "CG-DMR"], retained_col = "#ffca7b",
                 title_text = "CG-DMRs: Retained vs linked to SNP / SV / TIP"))
dev.off()
# =========================================================
# 14) association table per marker
# =========================================================
files <- list(
  SNP  = "/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/01.0_cis-trans_sig-associations.tsv",
  SV = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results/01.0_cis-trans_sig-associations.tsv",
  TIP = "/mnt/disk2/vibanez/otherAnalysis/11_DMR-GWAS-TIPs/ae_data-analysis/results/01.0_cis-trans_sig-associations.tsv"
)
# SNP-universe from your current analysis

norm_assoc <- function(dt, marker_set) {
  dt <- as.data.table(copy(dt))
  
  # standardize names
  nms <- names(dt)
  nms <- gsub("-", "_", nms)
  names(dt) <- nms
  
  # add missing columns if absent
  for (nm in c("chrSNP","start_loci","end_loci","meanBeta","meanSDbeta","minPvalue",
               "nSNPs","DMRpos","chrDMR","DMR","startDMR","lociType","lociPosition")) {
    if (!nm %in% names(dt)) dt[, (nm) := NA]
  }
  
  # force shared types
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
  
  dt[, marker_set := marker_set]
  
  dt[, .(chrSNP, start_loci, end_loci, meanBeta, meanSDbeta, minPvalue, nSNPs,
         DMRpos, chrDMR, DMR, startDMR, lociType, lociPosition, marker_set)]
}

assoc_list <- list(
  norm_assoc(fread(files[[1]]), "SNP"),
  norm_assoc(fread(files[[2]]), "SV"),
  norm_assoc(fread(files[[3]]), "TIP")
)

assoc_dt <- bind_rows(assoc_list)
#
toKeep<-status_long %>%
  filter(type=="sig")
status_long
assoc_dt <- assoc_dt %>%
  filter(DMRpos %in% dmr_universe$DMRid,DMRpos %in% toKeep$DMRid )
# define one QTL per row
mqtl_number <- assoc_dt %>%
  distinct(DMRid = DMRpos,
           dmr_class = DMR,
           marker_set,
           chrSNP,
           start_loci,
           end_loci) %>%
  count(DMRid, dmr_class, marker_set, name = "mQTL_number")

# classify number of QTLs
mqtl_full <- mqtl_number %>%
  mutate(class = case_when(
    mQTL_number == 0 ~ "0",
    mQTL_number == 1 ~ "1",
    mQTL_number <= 10 ~ "2-10",
    mQTL_number <= 20 ~ "11-20",
    mQTL_number <= 30 ~ "21-30",
    mQTL_number <= 40 ~ "31-40",
    mQTL_number <= 50 ~ "41-50",
    mQTL_number <= 100 ~ "51-100",
    TRUE ~ ">100"
  ))

mqtl_full$class <- factor(mqtl_full$class,levels = c("0", "1", "2-10", "11-20", "21-30", "31-40", "41-50", "51-100", ">100"))

# counts and proportions
plot_dt <- mqtl_full %>%
  count(dmr_class, marker_set, class, name = "count") %>%
  group_by(dmr_class, marker_set) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# plot
p <- ggplot(plot_dt, aes(x = proportion, y = class, fill = marker_set, color = marker_set)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7, width = 0.85) +
  facet_wrap(~ dmr_class, nrow = 1) +
  scale_fill_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  scale_color_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  labs(
    x = "Proportion of DMRs",
    y = "Number of QTLs",
    fill = "Marker set",
    color = "Marker set"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(OUT_QTL_PDF, p, width = 12, height = 8)

# total number of mQTLs per marker set
mqtl_total_marker <- mqtl_number %>%
  group_by(dmr_class, marker_set) %>%
  summarise(total_mQTL = sum(mQTL_number),
            n_DMR_with_mQTL = n_distinct(DMRid),
            .groups = "drop")

mqtl_total_marker$marker_set <- factor(mqtl_total_marker$marker_set,
  levels = c("SNP", "TIP", "SV"))

fwrite(as.data.table(mqtl_total_marker), OUT_QTL_TOTAL_TSV, sep = "\t")

p_total_mqtl <- ggplot(mqtl_total_marker,aes(x = marker_set, y = total_mQTL, fill = marker_set, color = marker_set)) +
  geom_col(alpha = 0.75, width = 0.7) +
  geom_text(aes(label = format(total_mQTL, big.mark = ",")), vjust = -0.35,size = 4  ) +
  facet_wrap(~ dmr_class, nrow = 1) +
  scale_fill_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  scale_color_manual(values = c(SNP = "#999999", TIP = "black", SV = "#977777")) +
  labs(x = "Marker set",y = "Total number of mQTLs",fill = "Marker set",color = "Marker set"  ) +
  theme_minimal() +
  theme(legend.position = "none",axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 11)  )

ggsave(OUT_QTL_TOTAL_PDF, p_total_mqtl, width = 10, height = 7)
