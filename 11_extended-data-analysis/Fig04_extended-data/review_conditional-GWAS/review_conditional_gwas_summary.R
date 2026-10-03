suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(eulerr)
})

# =========================
# CONFIG
# =========================
basedir_new <- "/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS/ad_outputs/"
basedir_old <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/"

outDir  <- "/mnt/disk2/vibanez/otherAnalysis/review_conditional-GWAS/"
outPlot <- outDir
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
dir.create(outPlot, recursive = TRUE, showWarnings = FALSE)

# Focus stratum
KEY_MM      <- "DMR"
KEY_KINSHIP <- "DMR-SNP"
KEY_INDEX   <- "aBN"
COVS <- c("SNP","TIP","SV")

# Threshold used for DMR-level significance
THR_DMR <- 1.312563e-07

# Colors
COL_COND <- c(base = "#c06a80", SNP = "#999999", TIP = "black", SV = "#977777")
COL_COV  <- c(SNP = "#999999", TIP = "black", SV = "#977777")
COL_LOST <- c(
  "Retained"         = "#820a86",
  "Lost_SNP"         = "#999999",
  "Lost_TIP"         = "black",
  "Lost_SV"          = "#977777",
  "Lost_SNP+TIP"     = "#666666",
  "Lost_SNP+SV"      = "#b08c8c",
  "Lost_TIP+SV"      = "#5f4b4b",
  "Lost_SNP+TIP+SV"  = "#404040"
)

# Output tables
f_base_tags         <- file.path(outDir, "review_cond_allDMR_base_tags.tsv")
f_trait_status      <- file.path(outDir, "review_cond_allDMR_trait_status.tsv")
f_dmr_long          <- file.path(outDir, "review_cond_allDMR_dmr_pvalues_long.tsv")
f_dmr_wide          <- file.path(outDir, "review_cond_allDMR_dmr_pvalues_wide.tsv")
f_summary           <- file.path(outDir, "review_cond_allDMR_summary.tsv")
f_lost_loci_wide    <- file.path(outDir, "review_cond_allDMR_lost_loci_wide.tsv")
f_lost_traits_wide  <- file.path(outDir, "review_cond_allDMR_lost_traits_wide.tsv")
f_partition_counts  <- file.path(outDir, "review_cond_allDMR_partition_counts.tsv")

# =========================
# HELPERS
# =========================
parse_tag <- function(stem) {
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  list(
    metabolite = parts[1],
    mm         = parts[2],
    kinship    = parts[3],
    index      = parts[4],
    covar      = ifelse(length(parts) >= 5, parts[5], "none")
  )
}

read_qtl <- function(path) {
  if (!file.exists(path)) return(data.table())
  dt <- tryCatch(fread(path, header = FALSE), error = function(e) data.table())
  if (nrow(dt) == 0) return(dt)
  setnames(dt, c("dmr_id","beta","beta_sd","p"))
  dt[p == 0, p := .Machine$double.xmin]
  dt
}

# Find only the result location / folder status
find_cond_result <- function(tag, covar) {
  stopifnot(covar %in% COVS)
  
  q_sig <- file.path(basedir_new, "sig",    sprintf("%s_%s.QTL",   tag, covar))
  q_non <- file.path(basedir_new, "nonSig", sprintf("%s_%s.QTL",   tag, covar))
  p_sig <- file.path(basedir_new, "sig",    sprintf("%s_%s.ps.gz", tag, covar))
  p_non <- file.path(basedir_new, "nonSig", sprintf("%s_%s.ps.gz", tag, covar))
  
  if (file.exists(q_sig) || file.exists(p_sig)) {
    return(list(
      folder = "sig",
      qtl    = if (file.exists(q_sig)) q_sig else NA_character_,
      ps     = if (file.exists(p_sig)) p_sig else NA_character_
    ))
  }
  
  if (file.exists(q_non) || file.exists(p_non)) {
    return(list(
      folder = "nonSig",
      qtl    = if (file.exists(q_non)) q_non else NA_character_,
      ps     = if (file.exists(p_non)) p_non else NA_character_
    ))
  }
  
  list(folder = NA_character_, qtl = NA_character_, ps = NA_character_)
}

extract_ps_many <- function(ps_gz, ids) {
  ids <- unique(ids[!is.na(ids)])
  if (length(ids) == 0) {
    return(data.table(
      dmr_id = character(),
      beta = numeric(),
      beta_sd = numeric(),
      p = numeric()
    ))
  }
  if (!file.exists(ps_gz)) {
    return(data.table(
      dmr_id = ids,
      beta = NA_real_,
      beta_sd = NA_real_,
      p = NA_real_
    ))
  }
  
  idfile <- tempfile()
  writeLines(ids, idfile)
  
  cmd <- sprintf(
    "zcat %s | awk 'NR==FNR{w[$1]=1; next} ($1 in w){print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}' %s -",
    shQuote(ps_gz), shQuote(idfile)
  )
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
  unlink(idfile)
  
  if (length(out) == 0) {
    return(data.table(
      dmr_id = ids,
      beta = NA_real_,
      beta_sd = NA_real_,
      p = NA_real_
    ))
  }
  
  dt <- fread(
    text = out,
    sep = "\t",
    header = FALSE,
    col.names = c("dmr_id","beta","beta_sd","p"),
    colClasses = c("character", "numeric", "numeric", "numeric")
  )
  
  dt[p == 0, p := .Machine$double.xmin]
  
  # ensure all requested ids are returned, even if missing in ps
  miss <- setdiff(ids, dt$dmr_id)
  if (length(miss) > 0) {
    dt <- rbind(
      dt,
      data.table(
        dmr_id = miss,
        beta = NA_real_,
        beta_sd = NA_real_,
        p = NA_real_
      ),
      fill = TRUE
    )
  }
  
  dt
}

make_combo_label3 <- function(snp, tip, sv, prefix = "Lost_") {
  out <- rep("Retained", length(snp))
  out[snp & !tip & !sv] <- paste0(prefix, "SNP")
  out[!snp & tip & !sv] <- paste0(prefix, "TIP")
  out[!snp & !tip & sv] <- paste0(prefix, "SV")
  out[snp & tip & !sv]  <- paste0(prefix, "SNP+TIP")
  out[snp & !tip & sv]  <- paste0(prefix, "SNP+SV")
  out[!snp & tip & sv]  <- paste0(prefix, "TIP+SV")
  out[snp & tip & sv]   <- paste0(prefix, "SNP+TIP+SV")
  out
}

combo_levels <- c(
  "Retained", "Lost_SNP", "Lost_TIP", "Lost_SV",
  "Lost_SNP+TIP", "Lost_SNP+SV", "Lost_TIP+SV", "Lost_SNP+TIP+SV"
)

combo_levels_traits <- c(
  "NONE", "SNP-only", "TIP-only", "SV-only",
  "SNP+TIP", "SNP+SV", "TIP+SV", "SNP+TIP+SV"
)

# =========================
# 1) BASELINE TAGS
# =========================
base_sig_dir <- file.path(basedir_old, "sig")
base_files   <- list.files(base_sig_dir, pattern = "\\.QTL$", full.names = FALSE)
base_stems   <- sub("\\.QTL$", "", base_files)

base_meta <- rbindlist(lapply(base_stems, function(st) {
  p <- parse_tag(st)
  data.table(
    tag        = st,
    metabolite = p$metabolite,
    mm         = p$mm,
    kinship    = p$kinship,
    index      = p$index
  )
}), fill = TRUE)

base_tags <- base_meta[mm == KEY_MM & kinship == KEY_KINSHIP & index == KEY_INDEX]

fwrite(base_tags, f_base_tags, sep = "\t")
message("Baseline traits (old sig): ", nrow(base_tags))

# =========================
# 2) TRAIT-LEVEL LOSS
# =========================
trait_status <- rbindlist(lapply(seq_len(nrow(base_tags)), function(i) {
  tag   <- base_tags$tag[i]
  metab <- base_tags$metabolite[i]
  
  rbindlist(lapply(COVS, function(cv) {
    info <- find_cond_result(tag, cv)
    data.table(
      metabolite      = metab,
      metabolite_code = sub("\\..*$", "", metab),
      tag             = tag,
      covar           = cv,
      cond_folder     = info$folder,
      trait_lost_all  = is.na(info$folder) || info$folder != "sig"
    )
  }))
}), fill = TRUE)

fwrite(trait_status, f_trait_status, sep = "\t")

lost_SNP <- sort(unique(trait_status[covar == "SNP" & trait_lost_all == TRUE, metabolite]))
lost_TIP <- sort(unique(trait_status[covar == "TIP" & trait_lost_all == TRUE, metabolite]))
lost_SV  <- sort(unique(trait_status[covar == "SV"  & trait_lost_all == TRUE, metabolite]))

cat("\n=== TRAITS THAT BECOME nonSig (folder-based) ===\n")
cat("Lost with SNP (n=", length(lost_SNP), "):\n", paste(lost_SNP, collapse = ", "), "\n\n", sep = "")
cat("Lost with TIP (n=", length(lost_TIP), "):\n", paste(lost_TIP, collapse = ", "), "\n\n", sep = "")
cat("Lost with SV  (n=", length(lost_SV ), "):\n", paste(lost_SV,  collapse = ", "), "\n\n", sep = "")

fit_all_traits <- euler(list(SNP = lost_SNP, TIP = lost_TIP, SV = lost_SV))
pdf(file.path(outPlot, "review_cond_01_trait_loss_venn_SNP_TIP_SV.pdf"), width = 5, height = 5)
plot(
  fit_all_traits,
  fills = c(COL_COV["SNP"], COL_COV["TIP"], COL_COV["SV"]),
  alpha = 0.35,
  edges = "black",
  quantities = TRUE,
  main = sprintf("Traits losing all DMR associations after conditioning (baseline n=%d)", nrow(base_tags))
)
dev.off()

# =========================
# 3) BASELINE DMR-QTL TABLE
# =========================
base_qtl <- rbindlist(lapply(seq_len(nrow(base_tags)), function(i) {
  tag   <- base_tags$tag[i]
  metab <- base_tags$metabolite[i]
  f     <- file.path(base_sig_dir, paste0(tag, ".QTL"))
  dt    <- read_qtl(f)
  if (nrow(dt) == 0) return(NULL)
  
  dt[, `:=`(
    tag          = tag,
    metabolite   = metab,
    p_base       = p,
    beta_base    = beta,
    beta_sd_base = beta_sd
  )]
  dt[, c("p","beta","beta_sd") := NULL]
  dt
}), fill = TRUE)

message("Baseline DMR-QTL rows: ", nrow(base_qtl))

# =========================
# 4) CONDITIONAL P-VALUES FROM ps.gz
# =========================
dmr_long <- rbindlist(lapply(COVS, function(cv) {
  x <- copy(base_qtl)
  x[, covar := cv]
  x
}), fill = TRUE)

dmr_long[, c("cond_folder", "ps_cond_path") := {
  info <- find_cond_result(tag[1], covar[1])
  list(info$folder, info$ps)
}, by = .(tag, covar)]

need_fill <- dmr_long[, .(tag, metabolite, covar, ps_cond_path, dmr_id)]

fills <- need_fill[, {
  ps <- ps_cond_path[1]
  if (is.na(ps) || !file.exists(ps)) {
    data.table(
      dmr_id       = dmr_id,
      beta_cond    = NA_real_,
      beta_sd_cond = NA_real_,
      p_cond       = NA_real_
    )
  } else {
    dt <- extract_ps_many(ps, dmr_id)
    setnames(dt, c("beta","beta_sd","p"), c("beta_cond","beta_sd_cond","p_cond"))
    dt
  }
}, by = .(tag, covar)]

setkey(dmr_long, tag, covar, dmr_id)
setkey(fills,    tag, covar, dmr_id)
dmr_long <- fills[dmr_long]

# final derived columns
dmr_long[, `:=`(
  base_sig  = (p_base <= THR_DMR),
  cond_sig  = (!is.na(p_cond) & p_cond <= THR_DMR),
  nlp_base  = -log10(p_base),
  nlp_cond  = ifelse(!is.na(p_cond), -log10(p_cond), NA_real_),
  dNL10p    = ifelse(!is.na(p_cond), (-log10(p_base)) - (-log10(p_cond)), NA_real_),
  lost_by_p = (p_base <= THR_DMR) & (is.na(p_cond) | p_cond > THR_DMR)
)]

# attach trait-level folder status
trait_status_dt <- unique(
  as.data.table(trait_status)[, .(metabolite, covar, trait_lost_all)],
  by = c("metabolite","covar")
)

setkey(trait_status_dt, metabolite, covar)
setkey(dmr_long, metabolite, covar)
dmr_long <- trait_status_dt[dmr_long]

dmr_long[, `:=`(
  covar          = as.character(covar),
  trait_lost_all = as.logical(trait_lost_all)
)]

setcolorder(dmr_long, c(
  "metabolite","covar","dmr_id","p_base","p_cond","beta_base","beta_sd_base",
  "beta_cond","beta_sd_cond","base_sig","cond_sig","lost_by_p",
  "nlp_base","nlp_cond","dNL10p","trait_lost_all","tag",
  "cond_folder","ps_cond_path"
))

fwrite(dmr_long, f_dmr_long, sep = "\t")

# =========================
# 5) DMR WIDE TABLE
# =========================
dmr_wide <- dcast( dmr_long,metabolite + dmr_id ~ covar,
  value.var = "p_cond", fun.aggregate = function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else min(x)})

base_p <- unique( as.data.table(base_qtl)[, .(metabolite, dmr_id, p_base)],
  by = c("metabolite","dmr_id"))

setkey(base_p, metabolite, dmr_id)
setkey(dmr_wide, metabolite, dmr_id)
dmr_wide <- base_p[dmr_wide]

for (cv in COVS) {
  if (!cv %in% names(dmr_wide)) dmr_wide[, (cv) := NA_real_]
}

setcolorder(dmr_wide, c("metabolite","dmr_id","p_base","SNP","TIP","SV"))
fwrite(dmr_wide, f_dmr_wide, sep = "\t")

# =========================
# 6) SUPPLEMENTARY TABLE JOIN
# =========================
supp_file <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/results/supplementary_tables-SuppTable_07.tsv"
if (file.exists(supp_file)) {
  supp <- fread(supp_file, sep = "\t", header = TRUE, skip = 1, fill = TRUE)
  pvals <- copy(dmr_wide)
  
  supp[, metabolite_key := paste0(`Metabolite name`, ".leaf-metabolites")]
  supp[, dmr_key := str_replace(`DMR name`, "^([^_]+)_([^_]+)_(.+)$", "\\1:\\3:\\2")]
  supp[, supp_row := .I]
  
  pvals[, metabolite_key := metabolite]
  pvals[, dmr_key := dmr_id]
  
  pvals_sub <- pvals[, .(
    metabolite_key, dmr_key,
    p_base,
    p_SNP = SNP,
    p_TIP = TIP,
    p_SV  = SV
  )]
  
  setkey(supp, metabolite_key, dmr_key)
  setkey(pvals_sub, metabolite_key, dmr_key)
  
  supp2 <- pvals_sub[supp]
  setorder(supp2, supp_row)
  supp2[, supp_row := NULL]
  
  fwrite(supp2, file.path(outDir, "supplementary_tables-SuppTable_07_with_pvalues.tsv"),sep = "\t", quote = FALSE)
  }

# =========================
# 7) SUMMARY TABLE
# =========================
summary_cov <- dmr_long[, .(
  n_traits_base          = uniqueN(metabolite),
  n_traits_lost          = uniqueN(metabolite[trait_lost_all == TRUE]),
  n_dmrs_base            = .N,
  n_dmrs_retained_sig    = sum(cond_sig, na.rm = TRUE),
  frac_dmrs_retained_sig = mean(cond_sig, na.rm = TRUE),
  pct_ge1                = mean(dNL10p >= 1, na.rm = TRUE),
  pct_ge2                = mean(dNL10p >= 2, na.rm = TRUE),
  median_dNL10p          = median(dNL10p, na.rm = TRUE)
), by = .(covar)]

summary_cov[, label := sprintf(
  "traits lost=%d/%d (%.1f%%); DMRs retained=%d/%d (%.1f%%); ≥10x=%.1f%%; ≥100x=%.1f%%; median dNL10p=%.2f",
  n_traits_lost, n_traits_base, 100 * n_traits_lost / n_traits_base,
  n_dmrs_retained_sig, n_dmrs_base, 100 * frac_dmrs_retained_sig,
  100 * pct_ge1, 100 * pct_ge2, median_dNL10p
)]

print(summary_cov)
fwrite(summary_cov, f_summary, sep = "\t")

cat("\n=== SUMMARY (DMR-level) ===\n")
print(summary_cov)

# =========================
# 8) PLOTS
# =========================
cap <- 20

dens_df <- rbindlist(list(
  dmr_long[, .(covar, nlp = pmin(nlp_base, cap), cond = "base")],
  dmr_long[, .(covar, nlp = pmin(nlp_cond, cap), cond = as.character(covar))]
), fill = TRUE)

dens_df[, cond  := factor(cond, levels = c("base","SNP","TIP","SV"))]
dens_df[, covar := factor(covar, levels = COVS)]

p_density <- ggplot(dens_df, aes(x = nlp, color = cond, fill = cond)) +
  geom_density(alpha = 0.20, linewidth = 1) +
  geom_vline(xintercept = -log10(THR_DMR), linetype = 2) +
  facet_wrap(~covar, nrow = 1) +
  scale_color_manual(values = COL_COND, limits = c("base","SNP","TIP","SV"), drop = FALSE) +
  scale_fill_manual(values = COL_COND, limits = c("base","SNP","TIP","SV"), drop = FALSE) +
  labs(
    x = expression(-log[10](p) ~ "(cap 20)"),
    y = "Density",
    title = "DMR-QTL p-values before vs after conditioning (baseline set)",
    subtitle = paste(summary_cov$covar, summary_cov$label, collapse = " | ")
  ) +
  theme_bw()

ggsave(file.path(outPlot, "review_cond_02_allDMR_density_pvalue_effect.pdf"),
       p_density, width = 32, height = 20, units = "cm")

p_violin <- ggplot( dmr_long[, .(covar, dNL10p)],
  aes(x = covar, y = pmin(pmax(dNL10p, -10), 10), fill = covar)) +
  geom_violin(trim = TRUE, alpha = 0.25) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = COL_COV, drop = FALSE) +
  labs(x = "Conditioning covariate", y = "dNL10p (cap ±10)",
    title = "Weakening of baseline DMR-QTL after conditioning (all baseline QTLs)"  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(outPlot, "review_cond_03_allDMR_violin_weakening.pdf"),
       p_violin, width = 24, height = 18, units = "cm")

p_ecdf <- ggplot(dmr_long, aes(x = dNL10p, color = covar)) +
  stat_ecdf(linewidth = 1) +
  coord_cartesian(xlim = c(-2, 5)) +
  geom_vline(xintercept = c(0, 1, 2), linetype = 2) +
  scale_color_manual(values = COL_COV, drop = FALSE) +
  labs(
    x = "dNL10p = (-log10 p_base) - (-log10 p_cond)",
    y = "ECDF",
    title = "Magnitude of weakening after conditioning (all baseline DMR-QTLs)",
    subtitle = "x=1 means ≥10× weaker; x=2 means ≥100× weaker."
  ) +
  theme_bw()

ggsave(file.path(outPlot, "review_cond_04_allDMR_ecdf_weakening.pdf"),
       p_ecdf, width = 24, height = 18, units = "cm")

p_scatter <- ggplot(dmr_long,aes(x = pmin(nlp_base, 20), y = pmin(nlp_cond, 20), color = cond_sig)) +
  geom_point(alpha = 0.35, size = 2.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~covar) +
  scale_color_manual(values = c(`TRUE` = "#3ac3c5", `FALSE` = "#ca3536")) +
  labs(
    x = expression(-log[10](p[base]) ~ "(cap 20)"),
    y = expression(-log[10](p[cond]) ~ "(cap 20)"),
    color = "DMR retained sig?",
    title = "Baseline DMR-QTL p-values before vs after conditioning (all baseline QTLs)"
  ) +
  theme_bw()

ggsave(file.path(outPlot, "review_cond_05_allDMR_scatter_pvalue_effect.pdf"),
       p_scatter, width = 28, height = 20, units = "cm")

# =========================
# 9) LOST LOCI OVERLAP
# =========================
dt <- copy(as.data.table(dmr_long))
dt[, covar := toupper(trimws(as.character(covar)))]
dt <- dt[covar %in% COVS]
dt[, locus := paste(metabolite, dmr_id, sep = "|")]

lost_tbl <- dt[, .(lost = any(lost_by_p, na.rm = TRUE)), by = .(covar, locus)]
lost_wide <- dcast(lost_tbl, locus ~ covar, value.var = "lost", fill = FALSE)

for (cv in COVS) {
  if (!cv %in% names(lost_wide)) lost_wide[, (cv) := FALSE]
}

lost_wide[, category := make_combo_label3(SNP, TIP, SV, prefix = "Lost_")]
lost_wide[, category := factor(category, levels = combo_levels)]
lost_wide[, c("metabolite","dmr_id") := tstrsplit(locus, "\\|", fixed = FALSE)]

fwrite(lost_wide, f_lost_loci_wide, sep = "\t")

plot_df <- lost_wide[, .N, by = category]
plot_df[, category := factor(category, levels = combo_levels)]

p_lost_bar <- ggplot(plot_df, aes(x = category, y = N, fill = category, color = category)) +
  geom_col() +
  geom_text(aes(label = N), vjust = -0.35, size = 3) +
  scale_fill_manual(values = COL_LOST, drop = FALSE) +
  scale_color_manual(values = COL_LOST, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "Number of loci",
    title = "Baseline DMR-QTL loci that lose significance after conditioning"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

ggsave(file.path(outPlot, "review_cond_06_lost_loci_overlap_barplot.pdf"),
       p_lost_bar, width = 24, height = 14, units = "cm")

fit_loci <- euler(list(
  SNP = lost_wide[SNP == TRUE, locus],
  TIP = lost_wide[TIP == TRUE, locus],
  SV  = lost_wide[SV  == TRUE, locus]
))
pdf(file.path(outPlot, "review_cond_06b_lost_loci_euler_SNP_TIP_SV.pdf"), width = 5, height = 5)
plot(
  fit_loci,
  fills = c(COL_COV["SNP"], COL_COV["TIP"], COL_COV["SV"]),
  alpha = 0.35,
  edges = "black",
  quantities = TRUE,
  main = "Lost DMR-QTL loci after conditioning"
)
dev.off()

# =========================
# 10) PARTITION ANALYSIS
# =========================
GENERAL_INFO <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/results/01.0_general-information.tsv"
if (file.exists(GENERAL_INFO)) {
  gi <- fread(GENERAL_INFO)
  
  dmr_base_traits <- unique(gi[
    mm == KEY_MM & kinship == KEY_KINSHIP & index == KEY_INDEX & resultType == "sig",
    metabolite
  ])
  
  snp_base_traits <- unique(gi[
    mm == "SNP" & kinship == "SNP" & index == KEY_INDEX & resultType == "sig",
    metabolite
  ])
  
  dmr_only_traits <- setdiff(dmr_base_traits, snp_base_traits)
  dmr_snp_traits  <- intersect(dmr_base_traits, snp_base_traits)
  
  ts <- as.data.table(trait_status)
  ts[, metabolite_code := ifelse(grepl("\\.leaf-metabolites$", metabolite),
                                 sub("\\..*$", "", metabolite), metabolite)]
  ts[, covar := toupper(trimws(as.character(covar)))]
  ts <- ts[covar %in% COVS]
  ts <- ts[metabolite_code %in% dmr_base_traits]
  
  lm <- ts[trait_lost_all == TRUE, .(metabolite = metabolite_code, covar)]
  
  if (nrow(lm) > 0) {
    lost_status_wide <- dcast(
      lm,
      metabolite ~ covar,
      value.var = "covar",
      fun.aggregate = function(x) length(x) > 0,
      fill = FALSE
    )
    
    for (cv in COVS) {
      if (!cv %in% names(lost_status_wide)) lost_status_wide[, (cv) := FALSE]
    }
    
    lost_status_wide[, partition := fifelse(
      metabolite %in% dmr_snp_traits, "DMR+SNP", "DMR-only"
    )]
    
    lost_status_wide[, lost_class := "NONE"]
    lost_status_wide[SNP & !TIP & !SV, lost_class := "SNP-only"]
    lost_status_wide[!SNP & TIP & !SV, lost_class := "TIP-only"]
    lost_status_wide[!SNP & !TIP & SV, lost_class := "SV-only"]
    lost_status_wide[SNP & TIP & !SV,  lost_class := "SNP+TIP"]
    lost_status_wide[SNP & !TIP & SV,  lost_class := "SNP+SV"]
    lost_status_wide[!SNP & TIP & SV,  lost_class := "TIP+SV"]
    lost_status_wide[SNP & TIP & SV,   lost_class := "SNP+TIP+SV"]
    
    lost_status_wide[, lost_class := factor(lost_class, levels = combo_levels_traits)]
    fwrite(lost_status_wide, f_lost_traits_wide, sep = "\t")
    
    counts_partition <- lost_status_wide[lost_class != "NONE",
                                         .(n_traits = .N), by = .(partition, lost_class)][order(partition, lost_class)]
    fwrite(counts_partition, f_partition_counts, sep = "\t")
    
    cat("\n=== LOST TRAITS BY PARTITION × CLASS ===\n")
    for (pt in unique(lost_status_wide$partition)) {
      cat("\n##", pt, "\n")
      for (cl in combo_levels_traits[combo_levels_traits != "NONE"]) {
        vv <- sort(lost_status_wide[partition == pt & lost_class == cl, metabolite])
        cat(cl, "(n=", length(vv), "): ", paste(vv, collapse = ", "), "\n", sep = "")
      }
    }
    
    venn_sets <- lapply(split(lost_status_wide, lost_status_wide$partition), function(dt0) {
      list(
        SNP = dt0[SNP == TRUE, metabolite],
        TIP = dt0[TIP == TRUE, metabolite],
        SV  = dt0[SV  == TRUE, metabolite]
      )
    })
    
    for (pt in names(venn_sets)) {
      sets <- venn_sets[[pt]]
      fit <- euler(sets)
      p<- plot(fit, fills = c(SNP = COL_COV["SNP"], TIP = COL_COV["TIP"], SV = COL_COV["SV"]),
               alpha = 0.35, edges = "black", quantities = TRUE, 
               main = sprintf("Lost traits (%s) after conditioning", pt))
      pdf(file.path(outPlot, sprintf("review_cond_07_trait_loss_venn_%s_SNP_TIP_SV.pdf",
                                     gsub("[^A-Za-z0-9]+", "_", pt))), width = 5, height = 5)
      print(p)
      dev.off()
    }
    
    p_partition <- ggplot(counts_partition, aes(x = lost_class, y = n_traits, fill = partition)) +
      geom_col(position = "dodge") +
      geom_text( aes(label = n_traits),position = position_dodge(width = 0.9), vjust = -0.3,size = 3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = NULL, y = "Number of traits", title = "Lost traits by partition and covariate combination") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    
    ggsave(file.path(outPlot, "review_cond_08_partition_lost_traits_barplot.pdf"), p_partition, width = 26, height = 14, units = "cm")
  }
}
#=== LOST TRAITS BY PARTITION × CLASS ===
  
## DMR+SNP 
# SNP-only(n=0): 
# TIP-only(n=0): 
# SV-only(n=0): 
# SNP+TIP(n=0): 
# SNP+SV(n=0): 
# TIP+SV(n=1): Tm003
# SNP+TIP+SV(n=2): Tm062, Tm221

## DMR-only 
# SNP-only(n=0): 
# TIP-only(n=1): Tm049
# SV-only(n=6): Tm065, Tm110, Tm129, Tm162, Tm205, Tm245
# SNP+TIP(n=1): Tm252
# SNP+SV(n=0): 
# TIP+SV(n=0): 
# SNP+TIP+SV(n=4): Tm193, Tm201, Tm203, Tm278
# =========================
# 11) OPTIONAL QUICK CONSOLE CHECKS
# =========================
cat("\nDone.\n")
cat("Saved:\n")
cat(" -", f_base_tags, "\n")
cat(" -", f_trait_status, "\n")
cat(" -", f_dmr_long, "\n")
cat(" -", f_dmr_wide, "\n")
cat(" -", f_summary, "\n")
cat(" -", f_lost_loci_wide, "\n")
cat(" -", f_lost_traits_wide, "\n")
cat(" -", f_partition_counts, "\n")
