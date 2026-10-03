#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(ggplot2)
})

# =========================================================
# 0) Paths and config
# =========================================================
outDir <- "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ad_heritability-merged"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
dmr_classes <- c("C-DMR", "CG-DMR")
chrs <- sprintf("ch%02d", 1:12)

# ---- 
cfg <- list(SNP = list(
    root = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results",
    layout = "flat",
    h2_types = c("sig", "nonSig","nonSig-GIF"),
    beta_types = c("sig")),
  TIP = list(root = "/mnt/disk2/vibanez/otherAnalysis/11_DMR-GWAS-TIPs/ab_results",
    layout = "nested",
    h2_types = c("sig", "nonSig","nonSig-GIF"),
    beta_types = c("sig")),
  SV = list(root = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results",
    layout = "nested",
    h2_types = c("sig", "nonSig", "nonSig-GIF"),
    beta_types = c("sig") )
)

h2_outfile      <- file.path(outDir, "00.1_merged-H2_all-markers.tsv")
beta_outfile    <- file.path(outDir, "00.2_merged-beta_sig_all-markers.tsv")
summary_outfile <- file.path(outDir, "00.3_H2_summary_all-markers.tsv")
plot_outfile    <- file.path(outDir, "00.4_H2_violin_all-markers.pdf")

ncore <- max(1L, min(12L, detectCores() - 1L))

# =========================================================
# 1) Helpers
# =========================================================
parse_dmr_info_from_basename <- function(fp) {
  bn <- basename(fp)
  id <- sub("\\..*$", "", bn)
  
  m <- regexec("^(ch\\d{2})_(C-DMR|CG-DMR)_(\\d+)$", id)
  rr <- regmatches(id, m)[[1]]
  
  if (length(rr) == 4) {
    list(
      DMRid = id,
      chr = rr[2],
      dmr_class = rr[3]
    )
  } else {
    list(
      DMRid = id,
      chr = NA_character_,
      dmr_class = NA_character_
    )
  }
}

read_reml_file <- function(fp, marker_set, dmr_class = NA_character_, chr = NA_character_, type = NA_character_) {
  df <- tryCatch(fread(fp, sep = "\t", fill = TRUE, header = FALSE, na.strings = "NA", nThread = 1),
    error = function(e) NULL)
  if (is.null(df)) return(NULL)
  
  vals <- unlist(df, use.names = FALSE)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NULL)
  
  info <- parse_dmr_info_from_basename(fp)
  if (is.na(dmr_class)) dmr_class <- info$dmr_class
  if (is.na(chr)) chr <- info$chr
  
  out <- as.data.table(as.list(vals))
  out[, DMRid := info$DMRid]
  out[, marker_set := marker_set]
  out[, dmr_class := dmr_class]
  out[, chr := chr]
  out[, type := type]
  
  # standardize first 6 columns if present
  nn <- ncol(out)
  std_names <- c("logVariance", "logLik", "delta", "sigmaGen", "sigmaError", "h2")
  n_original <- length(vals)
  k <- min(6, n_original)
  out
}

read_mqtl_file <- function(fp, marker_set, dmr_class = NA_character_, chr = NA_character_, type = NA_character_) {
  df <- tryCatch(
    fread(fp, sep = "\t", fill = TRUE, header = FALSE, na.strings = "NA", nThread = 1),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) < 4) return(NULL)
  
  info <- parse_dmr_info_from_basename(fp)
  if (is.na(dmr_class)) dmr_class <- info$dmr_class
  if (is.na(chr)) chr <- info$chr
  
  df <- df[, 1:4]
  setnames(df, c("marker_id", "beta", "betaSD", "pvalue"))
  df[, DMRid := info$DMRid]
  df[, marker_set := marker_set]
  df[, dmr_class := dmr_class]
  df[, chr := chr]
  df[, type := type]
  
  df
}

get_files_for_job <- function(marker_set, root, layout, dmr_class, chr, type, ext) {
  if (layout == "nested") {
    dir_i <- file.path(root, dmr_class, chr, type)
    if (!dir.exists(dir_i)) return(character())
    return(Sys.glob(file.path(dir_i, paste0("*.", ext))))
  }
  
  # flat layout, typically SNP
  dir_i <- file.path(root, type)
  if (!dir.exists(dir_i)) return(character())
  
  pat <- paste0("^", chr, "_", dmr_class, "_.*\\.", gsub("\\.", "\\\\.", ext), "$")
  fs <- list.files(dir_i, pattern = pat, full.names = TRUE)
  return(fs)
}

# =========================================================
# 2) Build jobs
# =========================================================
make_jobs <- function(cfg, which_types = c("h2", "beta")) {
  jobs <- list()
  
  for (marker_set in names(cfg)) {
    this <- cfg[[marker_set]]
    
    types <- if ("h2" %in% which_types) this$h2_types else this$beta_types
    
    jobs[[marker_set]] <- CJ(
      marker_set = marker_set,
      dmr_class = dmr_classes,
      chr = chrs,
      type = types,
      unique = TRUE
    )[, `:=`(
      root = this$root,
      layout = this$layout
    )]
  }
  
  rbindlist(jobs, use.names = TRUE, fill = TRUE)
}

jobs_h2 <- make_jobs(cfg, "h2")
jobs_beta <- make_jobs(cfg, "beta")

# =========================================================
# 3) H2 merge
# =========================================================
read_h2_job <- function(i, jobs) {
  marker_set <- jobs$marker_set[i]
  dmr_class  <- jobs$dmr_class[i]
  chr        <- jobs$chr[i]
  type       <- jobs$type[i]
  root       <- jobs$root[i]
  layout     <- jobs$layout[i]
  
  reml_files <- get_files_for_job(marker_set, root, layout, dmr_class, chr, type, "reml")
  if (length(reml_files) == 0) return(NULL)
  
  message("H2: ", marker_set, " / ", dmr_class, " / ", chr, " / ", type, " : ", length(reml_files), " files")
  
  out <- lapply(reml_files, function(fp) {
    read_reml_file(
      fp = fp,
      marker_set = marker_set,
      dmr_class = dmr_class,
      chr = chr,
      type = type
    )
  })
  
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) return(NULL)
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

h2_parts <- mclapply(
  seq_len(nrow(jobs_h2)),
  read_h2_job,
  jobs = jobs_h2,
  mc.cores = ncore
)

h2_parts <- Filter(Negate(is.null), h2_parts)
h2merged <- rbindlist(h2_parts, use.names = TRUE, fill = TRUE)

if ("h2" %in% names(h2merged)) {
  h2merged[, h2 := suppressWarnings(as.numeric(h2))]
}
if ("sigmaGen" %in% names(h2merged)) {
  h2merged[, sigmaGen := suppressWarnings(as.numeric(sigmaGen))]
}
if ("sigmaError" %in% names(h2merged)) {
  h2merged[, sigmaError := suppressWarnings(as.numeric(sigmaError))]
}

fwrite(h2merged, h2_outfile, sep = "\t", quote = FALSE)

# =========================================================
# 4) beta merge
# =========================================================
read_beta_job <- function(i, jobs) {
  marker_set <- jobs$marker_set[i]
  dmr_class  <- jobs$dmr_class[i]
  chr        <- jobs$chr[i]
  type       <- jobs$type[i]
  root       <- jobs$root[i]
  layout     <- jobs$layout[i]
  
  mqtl_files <- get_files_for_job(marker_set, root, layout, dmr_class, chr, type, "mQTL")
  if (length(mqtl_files) == 0) return(NULL)
  
  message("beta: ", marker_set, " / ", dmr_class, " / ", chr, " / ", type, " : ", length(mqtl_files), " files")
  
  out <- lapply(mqtl_files, function(fp) {
    read_mqtl_file(
      fp = fp,
      marker_set = marker_set,
      dmr_class = dmr_class,
      chr = chr,
      type = type
    )
  })
  
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) return(NULL)
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

beta_parts <- mclapply(
  seq_len(nrow(jobs_beta)),
  read_beta_job,
  jobs = jobs_beta,
  mc.cores = ncore
)

beta_parts <- Filter(Negate(is.null), beta_parts)
betaOut <- rbindlist(beta_parts, use.names = TRUE, fill = TRUE)

if ("beta" %in% names(betaOut)) {
  betaOut[, beta := suppressWarnings(as.numeric(beta))]
}
if ("betaSD" %in% names(betaOut)) {
  betaOut[, betaSD := suppressWarnings(as.numeric(betaSD))]
}
if ("pvalue" %in% names(betaOut)) {
  betaOut[, pvalue := suppressWarnings(as.numeric(pvalue))]
}

fwrite(betaOut, beta_outfile, sep = "\t", quote = FALSE)

# =========================================================
# 5) Summary table
# =========================================================
h2summary <- h2merged[
  !is.na(h2),
  .(
    n = .N,
    mean_h2 = mean(h2, na.rm = TRUE),
    median_h2 = median(h2, na.rm = TRUE),
    sd_h2 = sd(h2, na.rm = TRUE),
    q05 = quantile(h2, 0.05, na.rm = TRUE),
    q25 = quantile(h2, 0.25, na.rm = TRUE),
    q75 = quantile(h2, 0.75, na.rm = TRUE),
    q95 = quantile(h2, 0.95, na.rm = TRUE)
  ),
  by = .(marker_set, dmr_class, type)
][order(dmr_class, type, marker_set)]

fwrite(h2summary, summary_outfile, sep = "\t", quote = FALSE)

# =========================================================
# 6) Plot
# =========================================================
plot_dt <- copy(h2merged)[!is.na(h2)]
plot_dt[, marker_set := factor(marker_set, levels = c("SNP", "TIP", "SV"))]
plot_dt[, type := factor(type, levels = c("sig", "nonSig", "nonSig-GIF"))]
plot_dt[, dmr_class := factor(dmr_class, levels = c("C-DMR", "CG-DMR"))]

p <- ggplot(plot_dt, aes(x = marker_set, y = h2, fill = marker_set)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  facet_grid(dmr_class ~ type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(
    SNP = "#666666",
    TIP = "#6daedb",
    SV  = "#d95f5f"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Heritability (h2)", fill = "Marker set") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(plot_outfile, p, width = 28, height = 16, units = "cm")

message("Saved:")
message("  H2 merged   : ", h2_outfile)
message("  beta merged : ", beta_outfile)
message("  summary     : ", summary_outfile)
message("  plot        : ", plot_outfile)
