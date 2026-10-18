#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# =========================================================
# Compare raw BIseq and ONT methylation on the same BIseq DMRs
# Output:
#   1) joined methylation table: BIseq_* and ONT_* columns
#   2) BIseq-vs-ONT same-sample correlations
#   3) full sample corplot with BIseq and ONT together
# =========================================================

# -----------------------------
# Paths
# -----------------------------
ONT_data_dir <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ai_DMR-classification/bd_merged-methylation/"
ONT_data_loci <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ai_DMR-classification/bb_get-collapsed-loci/"

Biseq_data_dir <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation/cervil-denovo/"
Biseq_data_loci <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree/cervil-denovo/"

# BIseq DMRs lifted from SL2.5 to polished/ONT coordinates
Biseq_sl25_denovo_dir <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation/cervil-sl25/"

outDir <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/al_DMR-gen-var_denovo-pedigree_comparison/BIseq_ONT_raw_correlation/"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Parameters
# -----------------------------
samples <- c("G0-C3", "G0-C6", "G0-C9", "G5-C1", "G5-C2", "G5-C3")
dmr_types <- c("C-DMR", "CG-DMR")
method_cor <- "spearman"

# =========================================================
# Main loop
# =========================================================
all_joined <- list()
all_same_sample_cor <- list()

for (dmr_type in dmr_types) {

  message("\n==============================")
  message("Processing ", dmr_type)
  message("==============================")

  # 1) Read lifted BIseq DMRs in polished/ONT coordinates
    lifted_file <- file.path(Biseq_sl25_denovo_dir, paste0(dmr_type, ".SL25_to_polished.bed"))
  if (!file.exists(lifted_file)) stop("Missing lifted file: ", lifted_file)

  lifted <- fread(lifted_file, header = FALSE, sep = "\t", fill = TRUE)
  if (ncol(lifted) < 4) stop("Lifted file has fewer than 4 columns: ", lifted_file)

  lifted <- lifted[, .(
    DMR = dmr_type,
    chr = as.character(V1),
    start = as.integer(V2),
    end = as.integer(V3),
    dmr_id = as.character(V4)
  )]
  lifted <- unique(lifted)
  lifted[, lifted_row_id := .I]

  fwrite(
    lifted,
    file = file.path(outDir, paste0("00_", dmr_type, "_lifted_BIseq_DMRs_used_for_matching.tsv")),
    sep = "\t",
    quote = FALSE
  )

  message("Lifted BIseq DMRs: ", nrow(lifted))

    # 2) Read ONT methylation matrix
  ont_files <- list.files(
    ONT_data_dir,
    pattern = paste0("^", dmr_type, "_all.*\\.bed(\\.gz)?$"),
    full.names = TRUE
  )
  if (length(ont_files) == 0) stop("No ONT methylation BED found for: ", dmr_type)
  ont_file <- sort(ont_files)[1]

  ont_header_file <- file.path(ONT_data_dir, paste0(dmr_type, "_all.colNames.tsv"))
  if (!file.exists(ont_header_file)) stop("Missing ONT header: ", ont_header_file)

  ont_header <- fread(ont_header_file, header = FALSE)$V1
  ont <- fread(ont_file, header = FALSE, sep = "\t", fill = TRUE, na.strings = c("", "NA"))

  ont_names <- c("chr", "start", "end", ont_header)
  if (length(ont_names) != ncol(ont)) {
    stop("ONT column mismatch for ", ont_file,
         "\nExpected ", length(ont_names), " columns; observed ", ncol(ont))
  }
  setnames(ont, ont_names)

  ont_samples <- intersect(samples, colnames(ont))
  if (length(ont_samples) == 0) stop("No selected samples found in ONT file: ", ont_file)

  ont <- ont[, c("chr", "start", "end", ont_samples), with = FALSE]
  ont[, `:=`(
    DMR = dmr_type,
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end)
  )]
  ont[, (ont_samples) := lapply(.SD, as.numeric), .SDcols = ont_samples]

  # Optional: keep collapsed ONT loci, as in your previous scripts
  ont_loci_file <- file.path(ONT_data_loci, paste0(dmr_type, "-loci_collapsed.bed"))
  if (file.exists(ont_loci_file)) {
    ont_loci <- fread(ont_loci_file, header = FALSE, sep = "\t", fill = TRUE)
    setnames(ont_loci, c("V1", "V2"), c("chr", "start"))
    ont_loci <- unique(ont_loci[, .(
      chr = as.character(chr),
      start = as.integer(start)
    )])
    ont <- merge(ont, ont_loci, by = c("chr", "start"), all = FALSE)
  }

  ont[, ont_locus_id := paste(chr, start, end, sep = "_")]
  ont[, ont_row_id := .I]

  fwrite(
    head(ont, 100),
    file = file.path(outDir, paste0("01_", dmr_type, "_ONT_loaded_head100.tsv")),
    sep = "\t",
    quote = FALSE
  )

  message("ONT methylation rows loaded: ", nrow(ont))

  
  # 3) Read BIseq methylation matrix
  
  biseq_files <- list.files(
    Biseq_data_dir,
    pattern = paste0("^", dmr_type, "_all.*\\.bed(\\.gz)?$"),
    full.names = TRUE
  )
  if (length(biseq_files) == 0) stop("No BIseq methylation BED found for: ", dmr_type)
  biseq_file <- sort(biseq_files)[1]

  biseq_header_file <- file.path(Biseq_data_dir, paste0(dmr_type, "_all.colNames.tsv"))
  if (!file.exists(biseq_header_file)) stop("Missing BIseq header: ", biseq_header_file)

  biseq_header <- fread(biseq_header_file, header = FALSE)$V1
  biseq <- fread(biseq_file, header = FALSE, sep = "\t", fill = TRUE, na.strings = c("", "NA"))

  biseq_names <- c("chr", "start", "end", biseq_header)
  if (length(biseq_names) != ncol(biseq)) {
    stop("BIseq column mismatch for ", biseq_file,
         "\nExpected ", length(biseq_names), " columns; observed ", ncol(biseq))
  }
  setnames(biseq, biseq_names)

  biseq_samples <- intersect(samples, colnames(biseq))
  if (length(biseq_samples) == 0) stop("No selected samples found in BIseq file: ", biseq_file)

  biseq <- biseq[, c("chr", "start", "end", biseq_samples), with = FALSE]
  biseq[, `:=`(
    DMR = dmr_type,
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end)
  )]
  biseq[, (biseq_samples) := lapply(.SD, as.numeric), .SDcols = biseq_samples]

  # Optional: keep collapsed BIseq loci, as in your previous scripts
  biseq_loci_file <- file.path(Biseq_data_loci, paste0(dmr_type, "-loci_collapsed.bed"))
  if (file.exists(biseq_loci_file)) {
    biseq_loci <- fread(biseq_loci_file, header = FALSE, sep = "\t", fill = TRUE)
    setnames(biseq_loci, c("V1", "V2"), c("chr", "start"))
    biseq_loci <- unique(biseq_loci[, .(
      chr = as.character(chr),
      start = as.integer(start)
    )])
    biseq <- merge(biseq, biseq_loci, by = c("chr", "start"), all = FALSE)
  }

  biseq[, biseq_locus_id := paste(chr, start, end, sep = "_")]
  biseq[, biseq_row_id := .I]

  fwrite(
    head(biseq, 100),
    file = file.path(outDir, paste0("01_", dmr_type, "_BIseq_loaded_head100.tsv")),
    sep = "\t",
    quote = FALSE
  )

  message("BIseq methylation rows loaded: ", nrow(biseq))

  
  # 4) Match ONT methylation to lifted BIseq DMRs
  
  lifted_for_ont <- copy(lifted)
  setkey(lifted_for_ont, DMR, chr, start, end)
  setkey(ont, DMR, chr, start, end)

  ont_hits <- foverlaps(
    ont,
    lifted_for_ont,
    by.x = c("DMR", "chr", "start", "end"),
    by.y = c("DMR", "chr", "start", "end"),
    type = "any",
    nomatch = 0L
  )

  message("ONT overlap rows with lifted BIseq DMRs: ", nrow(ont_hits))

  if (nrow(ont_hits) == 0) {
    warning("No ONT overlap found for ", dmr_type)
    next
  }

  ont_ag <- ont_hits[, c(
    list(
      n_ont_loci = .N,
      ont_loci = paste(unique(ont_locus_id), collapse = ";")
    ),
    lapply(.SD, function(z) {
      z <- as.numeric(z)
      if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
    })
  ), by = .(DMR, chr, start, end, dmr_id, lifted_row_id), .SDcols = ont_samples]

  setnames(ont_ag, ont_samples, paste0("ONT_", ont_samples))

  fwrite(
    ont_ag,
    file = file.path(outDir, paste0("02_", dmr_type, "_ONT_aggregated_on_BIseq_DMRs.tsv")),
    sep = "\t",
    quote = FALSE
  )

  
  # 5) Match BIseq methylation to the same lifted BIseq DMRs
  
  lifted_for_biseq <- copy(lifted)
  setkey(lifted_for_biseq, DMR, chr, start, end)
  setkey(biseq, DMR, chr, start, end)

  biseq_hits <- foverlaps(
    biseq,
    lifted_for_biseq,
    by.x = c("DMR", "chr", "start", "end"),
    by.y = c("DMR", "chr", "start", "end"),
    type = "any",
    nomatch = 0L
  )

  message("BIseq overlap rows with lifted BIseq DMRs: ", nrow(biseq_hits))

  if (nrow(biseq_hits) == 0) {
    warning("No BIseq overlap found for ", dmr_type,
            ". Check whether Biseq_data_dir and lifted file are in the same coordinates.")
    next
  }

  biseq_ag <- biseq_hits[, c(
    list(
      n_biseq_loci = .N,
      biseq_loci = paste(unique(biseq_locus_id), collapse = ";")
    ),
    lapply(.SD, function(z) {
      z <- as.numeric(z)
      if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
    })
  ), by = .(DMR, chr, start, end, dmr_id, lifted_row_id), .SDcols = biseq_samples]

  setnames(biseq_ag, biseq_samples, paste0("BIseq_", biseq_samples))

  fwrite(
    biseq_ag,
    file = file.path(outDir, paste0("02_", dmr_type, "_BIseq_aggregated_on_BIseq_DMRs.tsv")),
    sep = "\t",
    quote = FALSE
  )

  
  # 6) Join BIseq and ONT into one visible comparison table
  
  joined <- merge(
    biseq_ag,
    ont_ag,
    by = c("DMR", "chr", "start", "end", "dmr_id", "lifted_row_id"),
    all = FALSE
  )

  setorder(joined, DMR, chr, start, end)
  message("Joined BIseq + ONT DMRs: ", nrow(joined))

  fwrite(
    joined,
    file = file.path(outDir, paste0("03_", dmr_type, "_BIseq_ONT_joined_methylation_table.tsv")),
    sep = "\t",
    quote = FALSE
  )

  # Print a small visible table in the terminal
  message("\nHead of joined table for ", dmr_type, ":")
  print(head(joined, 10))

  all_joined[[dmr_type]] <- joined

  
  # 7) Same-sample BIseq-vs-ONT correlations
  
  same_cor <- data.table()

  for (s in samples) {
    bcol <- paste0("BIseq_", s)
    ocol <- paste0("ONT_", s)

    if (!all(c(bcol, ocol) %in% colnames(joined))) next

    ok <- complete.cases(joined[[bcol]], joined[[ocol]])

    same_cor <- rbind(
      same_cor,
      data.table(
        DMR = dmr_type,
        sample = s,
        n = sum(ok),
        spearman = ifelse(sum(ok) >= 3, cor(joined[[bcol]][ok], joined[[ocol]][ok], method = "spearman"), NA_real_),
        pearson = ifelse(sum(ok) >= 3, cor(joined[[bcol]][ok], joined[[ocol]][ok], method = "pearson"), NA_real_)
      )
    )
  }

  fwrite(
    same_cor,
    file = file.path(outDir, paste0("04_", dmr_type, "_same_sample_BIseq_vs_ONT_correlations.tsv")),
    sep = "\t",
    quote = FALSE
  )

  all_same_sample_cor[[dmr_type]] <- same_cor

  message("\nSame-sample BIseq vs ONT correlations for ", dmr_type, ":")
  print(same_cor)

  # Scatter plot: BIseq vs ONT for each sample
  plot_long <- rbindlist(lapply(samples, function(s) {
    bcol <- paste0("BIseq_", s)
    ocol <- paste0("ONT_", s)
    if (!all(c(bcol, ocol) %in% colnames(joined))) return(NULL)
    data.table(
      DMR = dmr_type,
      sample = s,
      BIseq = joined[[bcol]],
      ONT = joined[[ocol]]
    )
  }), fill = TRUE)

  plot_long <- plot_long[!is.na(BIseq) & !is.na(ONT)]

  if (nrow(plot_long) > 0) {
    p_scatter <- ggplot(plot_long, aes(x = BIseq, y = ONT)) +
      geom_point(alpha = 0.25, size = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      facet_wrap(~sample, ncol = 3) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_minimal() +
      labs(
        title = paste0(dmr_type, ": raw methylation BIseq vs ONT"),
        x = "BIseq methylation",
        y = "ONT methylation"
      )

    ggsave(
      file.path(outDir, paste0("04_", dmr_type, "_BIseq_vs_ONT_scatter_by_sample.pdf")),
      p_scatter,
      width = 9,
      height = 6
    )
  }

  
  # 8) Full corplot: BIseq_* and ONT_* together
  
  biseq_cols <- paste0("BIseq_", samples)
  ont_cols <- paste0("ONT_", samples)
  cor_cols <- c(biseq_cols, ont_cols)
  cor_cols <- cor_cols[cor_cols %in% colnames(joined)]

  mat <- as.matrix(joined[, ..cor_cols])
  rownames(mat) <- make.unique(joined$dmr_id)

  # Keep rows with enough information for pairwise correlations
  mat <- mat[rowSums(!is.na(mat)) >= 4, , drop = FALSE]

  if (nrow(mat) >= 3) {
    cor_mat <- cor(mat, use = "pairwise.complete.obs", method = method_cor)

    cor_dt <- as.data.table(as.table(cor_mat))
    setnames(cor_dt, c("sample1", "sample2", "correlation"))

    fwrite(
      cor_dt,
      file = file.path(outDir, paste0("05_", dmr_type, "_BIseq_ONT_full_sample_correlation_", method_cor, ".tsv")),
      sep = "\t",
      quote = FALSE
    )

    pdf(
      file.path(outDir, paste0("05_", dmr_type, "_BIseq_ONT_full_sample_correlation_", method_cor, ".pdf")),
      width = 8,
      height = 7
    )

    pheatmap(
      cor_mat,
      color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      display_numbers = TRUE,
      number_format = "%.2f",
      main = paste0(dmr_type, ": BIseq + ONT sample correlation"),
      border_color = NA
    )

    dev.off()
  } else {
    warning("Not enough rows for full correlation plot: ", dmr_type)
  }
}


# 9) Combined outputs for both DMR classes
# =========================================================
if (length(all_joined) > 0) {
  joined_all <- rbindlist(all_joined, fill = TRUE)

  fwrite(
    joined_all,
    file = file.path(outDir, "03_ALL_BIseq_ONT_joined_methylation_table.tsv"),
    sep = "\t",
    quote = FALSE
  )

  message("\nCombined joined table:")
  print(joined_all[, .N, by = DMR])
  print(head(joined_all, 10))

  biseq_cols <- paste0("BIseq_", samples)
  ont_cols <- paste0("ONT_", samples)
  cor_cols <- c(biseq_cols, ont_cols)
  cor_cols <- cor_cols[cor_cols %in% colnames(joined_all)]

  mat_all <- as.matrix(joined_all[, ..cor_cols])
  rownames(mat_all) <- make.unique(paste(joined_all$DMR, joined_all$dmr_id, sep = "_"))
  mat_all <- mat_all[rowSums(!is.na(mat_all)) >= 4, , drop = FALSE]

  if (nrow(mat_all) >= 3) {
    cor_all <- cor(mat_all, use = "pairwise.complete.obs", method = method_cor)

    cor_all_dt <- as.data.table(as.table(cor_all))
    setnames(cor_all_dt, c("sample1", "sample2", "correlation"))

    fwrite(
      cor_all_dt,
      file = file.path(outDir, paste0("05_ALL_BIseq_ONT_full_sample_correlation_", method_cor, ".tsv")),
      sep = "\t",
      quote = FALSE
    )

    pdf(
      file.path(outDir, paste0("05_ALL_BIseq_ONT_full_sample_correlation_", method_cor, ".pdf")),
      width = 8,
      height = 7
    )

    pheatmap(
      cor_all,
      color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      display_numbers = TRUE,
      number_format = "%.2f",
      main = "ALL DMRs: BIseq + ONT sample correlation",
      border_color = NA
    )

    dev.off()
  }
}

if (length(all_same_sample_cor) > 0) {
  same_cor_all <- rbindlist(all_same_sample_cor, fill = TRUE)

  fwrite(
    same_cor_all,
    file = file.path(outDir, "04_ALL_same_sample_BIseq_vs_ONT_correlations.tsv"),
    sep = "\t",
    quote = FALSE
  )

  message("\nAll same-sample BIseq vs ONT correlations:")
  print(same_cor_all)
}

message("\nDone. Outputs written to: ", outDir)
