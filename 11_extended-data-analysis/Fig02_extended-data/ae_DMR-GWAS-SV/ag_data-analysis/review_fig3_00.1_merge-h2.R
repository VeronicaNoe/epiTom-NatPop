suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

wd <- "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results"
outDir <- "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results/"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

dmr_classes <- c("CG-DMR", "C-DMR")
chrs <- sprintf("ch%02d", 1:12)
types_h2 <- c("sig", "nonSig", "nonSig-GIF")
types_beta <- c("sig")

h2_outfile   <- file.path(outDir, "00.1_merged-H2_allDMRs.tsv")
beta_outfile <- file.path(outDir, "00.1_merged-beta-values_sig-DMRs.tsv")

jobs_h2 <- CJ(dmr_class = dmr_classes, chr = chrs, type = types_h2, unique = TRUE)
jobs_beta <- CJ(dmr_class = dmr_classes, chr = chrs, type = types_beta, unique = TRUE)

ncore <- min(8L, detectCores() - 1L)

read_h2_job <- function(d, ch, tp) {
  dir_i <- file.path(wd, d, ch, tp)
  if (!dir.exists(dir_i)) return(NULL)
  
  reml_files <- Sys.glob(file.path(dir_i, "*.reml"))
  if (length(reml_files) == 0) return(NULL)
  
  message("H2: ", d, " / ", ch, " / ", tp, " : ", length(reml_files), " files")
  
  h2_list <- lapply(reml_files, function(fp) {
    df <- tryCatch(
      fread(fp, sep = "\t", fill = TRUE, header = FALSE, na.strings = "NA", nThread = 1),
      error = function(e) NULL
    )
    if (is.null(df)) return(NULL)
    
    vals <- unlist(df, use.names = FALSE)
    vals <- vals[!is.na(vals)]
    
    out <- as.data.table(as.list(vals))
    out[, DMR := sub("\\..*$", "", basename(fp))]
    out[, dmr_class := d]
    out[, chr := ch]
    out[, type := tp]
    out
  })
  
  h2_list <- Filter(Negate(is.null), h2_list)
  if (length(h2_list) == 0) return(NULL)
  
  res <- rbindlist(h2_list, use.names = TRUE, fill = TRUE)
  
  if (ncol(res) >= 10) {
    setnames(
      res,
      old = names(res)[1:6],
      new = c("logVariance", "log", "delta", "sigmaGen", "sigmaError", "h2")
    )
  }
  
  res
}

h2_parts <- mclapply(
  seq_len(nrow(jobs_h2)),
  function(i) read_h2_job(jobs_h2$dmr_class[i], jobs_h2$chr[i], jobs_h2$type[i]),
  mc.cores = ncore
)

h2_parts <- Filter(Negate(is.null), h2_parts)
h2merged <- rbindlist(h2_parts, use.names = TRUE, fill = TRUE)

fwrite(h2merged, h2_outfile, sep = "\t", quote = FALSE)

read_beta_job <- function(d, ch, tp) {
  dir_i <- file.path(wd, d, ch, tp)
  if (!dir.exists(dir_i)) return(NULL)
  
  mqtl_files <- Sys.glob(file.path(dir_i, "*.mQTL"))
  if (length(mqtl_files) == 0) return(NULL)
  
  message("beta: ", d, " / ", ch, " / ", tp, " : ", length(mqtl_files), " files")
  
  beta_list <- lapply(mqtl_files, function(fp) {
    df <- tryCatch(
      fread(fp, sep = "\t", fill = TRUE, header = FALSE, na.strings = "NA", nThread = 1),
      error = function(e) NULL
    )
    if (is.null(df) || ncol(df) < 4) return(NULL)
    
    df <- df[, 1:4]
    setnames(df, c("SNPid", "beta", "betaSD", "pvalues"))
    df[, DMRid := sub("\\..*$", "", basename(fp))]
    df[, dmr_class := d]
    df[, chr := ch]
    df[, type := tp]
    df
  })
  
  beta_list <- Filter(Negate(is.null), beta_list)
  if (length(beta_list) == 0) return(NULL)
  
  rbindlist(beta_list, use.names = TRUE, fill = TRUE)
}

beta_parts <- mclapply(
  seq_len(nrow(jobs_beta)),
  function(i) read_beta_job(jobs_beta$dmr_class[i], jobs_beta$chr[i], jobs_beta$type[i]),
  mc.cores = ncore
)

beta_parts <- Filter(Negate(is.null), beta_parts)
betaOut <- rbindlist(beta_parts, use.names = TRUE, fill = TRUE)

fwrite(betaOut, beta_outfile, sep = "\t", quote = FALSE)
