# =========================================================
# Plot ChIP-seq and ATAC-seq continuous signal over DMRs
# grouped by fUM bin and DMR annotation
#
# Input:
#   - DMR BED4 files
#   - fUM frequency files
#   - ChIP bigWig signal summarized by multiBigwigSummary
#   - ATAC bigWig signal summarized by multiBigwigSummary
#   - region annotation BED4 files
#
# Output:
#   - median ChIP signal by fUM bin / region / mark
#   - median ATAC signal by fUM bin / region
#   - combined plots with ChIP and ATAC signal
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================================================
# Paths
# =========================================================

OUTROOT <- "/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq"

outPlot <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/plots/"
outDir  <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/results/"

dir.create(outPlot, showWarnings = FALSE, recursive = TRUE)
dir.create(outDir,  showWarnings = FALSE, recursive = TRUE)

DMRSETS <- c("CG-DMR", "C-DMR")

# =========================================================
# Helper functions
# =========================================================

read_ids <- function(f) {
  if (!file.exists(f)) return(character())
  x <- fread(f, header = FALSE)
  if (nrow(x) == 0) return(character())
  as.character(x[[4]])
}

read_deeptools_signal <- function(signal_file, bed4_file, signal_names) {

  if (!file.exists(signal_file)) {
    stop("Missing signal file: ", signal_file)
  }

  if (!file.exists(bed4_file)) {
    stop("Missing BED4 file: ", bed4_file)
  }

  X <- fread(signal_file)

  # deepTools outRawCounts usually has:
  # chr start end sample1 sample2 ...
  # Sometimes the first columns may be named differently.
  # We force the expected names.
  n_expected <- 3 + length(signal_names)

  if (ncol(X) != n_expected) {
    stop(
      "Unexpected number of columns in ", signal_file,
      ". Expected ", n_expected, " but found ", ncol(X), "."
    )
  }

  setnames(X, c("chr", "start", "end", signal_names))

  X[, start := as.integer(start)]
  X[, end   := as.integer(end)]

  for (cc in signal_names) {
    X[, (cc) := as.numeric(get(cc))]
  }

  B <- fread(bed4_file, header = FALSE)
  setnames(B, c("chr", "start", "end", "id"))

  B[, start := as.integer(start)]
  B[, end   := as.integer(end)]
  B[, id    := as.character(id)]

  X_id <- merge(
    X,
    B,
    by = c("chr", "start", "end"),
    all.x = TRUE
  )

  if (any(is.na(X_id$id))) {
    warning(
      "Some signal rows could not be matched to DMR IDs in: ",
      basename(signal_file)
    )
  }

  X_id
}

rescale_to_range <- function(x, new_min = 0, new_max = 1) {
  old_min <- min(x, na.rm = TRUE)
  old_max <- max(x, na.rm = TRUE)

  if (!is.finite(old_min) || !is.finite(old_max) || old_min == old_max) {
    return(rep(NA_real_, length(x)))
  }

  new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min)
}

# =========================================================
# Main loop
# =========================================================

for (dmrset in DMRSETS) {

  message("Processing ", dmrset)

  # -------------------------------------------------------
  # Files
  # -------------------------------------------------------

  freq_file <- sprintf(
    "%s/ad_UM-frequency/master/%s.methfreq.tsv.gz",
    OUTROOT, dmrset
  )

  bed4_file <- file.path(
    OUTROOT,
    "aa_prepare-public-data/master",
    paste0(dmrset, ".M82.bed4")
  )

  chip_signal_file <- file.path(
    OUTROOT,
    "aa_prepare-public-data/signal_by_DMR",
    paste0(dmrset, ".ChIP_signal_by_DMR.tsv")
  )

  atac_signal_file <- file.path(
    OUTROOT,
    "aa_prepare-public-data/signal_by_DMR",
    paste0(dmrset, ".ATAC_signal_by_DMR.tsv")
  )

  # -------------------------------------------------------
  # Read DMR coordinates and fUM
  # -------------------------------------------------------

  B <- fread(bed4_file, header = FALSE)
  setnames(B, c("chr", "start", "end", "id"))

  B[, id := as.character(id)]
  B[, start := as.integer(start)]
  B[, end   := as.integer(end)]

  F <- fread(freq_file)
  setnames(
    F,
    c("id", "chr_freq", "pos_freq", "end_freq", "n_um", "n_m", "n_obs", "fUM", "fMeth")
  )

  F[, id := as.character(id)]

  D <- merge(
    B,
    F[, .(id, n_obs, fUM, fMeth)],
    by = "id",
    all.x = TRUE
  )

  # -------------------------------------------------------
  # Add region annotations
  # -------------------------------------------------------

  REGDIR <- file.path(OUTROOT, "ac_annotate-DMRs-M82/dmrs_by_region")

  f_prom   <- file.path(REGDIR, sprintf("%s.promoter.bed4", dmrset))
  f_gene   <- file.path(REGDIR, sprintf("%s.gene.bed4", dmrset))
  f_geneTE <- file.path(REGDIR, sprintf("%s.gene-TE.bed4", dmrset))
  f_te     <- file.path(REGDIR, sprintf("%s.TE.bed4", dmrset))
  f_inter  <- file.path(REGDIR, sprintf("%s.intergenic.bed4", dmrset))

  ids_prom   <- read_ids(f_prom)
  ids_gene   <- read_ids(f_gene)
  ids_geneTE <- read_ids(f_geneTE)
  ids_te     <- read_ids(f_te)
  ids_inter  <- read_ids(f_inter)

  ids_gene_noTE <- setdiff(ids_gene, ids_geneTE)

  D[, dmr_loc := NA_character_]

  D[id %in% ids_prom, dmr_loc := "promoter"]
  D[is.na(dmr_loc) & id %in% ids_geneTE,    dmr_loc := "gene-TE"]
  D[is.na(dmr_loc) & id %in% ids_gene_noTE, dmr_loc := "gene"]
  D[is.na(dmr_loc) & id %in% ids_te,        dmr_loc := "TE"]
  D[is.na(dmr_loc) & id %in% ids_inter,     dmr_loc := "intergenic"]

  D[, dmr_loc := factor(
    dmr_loc,
    levels = c("promoter", "gene-TE", "gene", "TE", "intergenic")
  )]

  print(D[, .N, by = dmr_loc][order(-N)])

  # -------------------------------------------------------
  # fUM bins
  # -------------------------------------------------------

  breaks <- seq(0, 1, by = 0.1)

  labels <- paste0(
    sprintf("%.1f", breaks[-length(breaks)]),
    "-",
    sprintf("%.1f", breaks[-1])
  )

  D[, fUM_bin := cut(
    fUM,
    breaks = breaks,
    include.lowest = TRUE,
    right = TRUE,
    labels = labels
  )]

  D[, fUM_bin := factor(fUM_bin, levels = labels)]

  # -------------------------------------------------------
  # Read ChIP and ATAC signal tables
  # -------------------------------------------------------

  chip_names <- c(
    "H3K4me3", "H3K9ac", "H3K27ac", "H3K36me3",
    "H3K9me2", "Pol2", "H3K27me3"
  )

  atac_names <- c("ATAC_rep1", "ATAC_rep2")

  CHIP <- read_deeptools_signal(
    signal_file = chip_signal_file,
    bed4_file = bed4_file,
    signal_names = chip_names
  )

  ATAC <- read_deeptools_signal(
    signal_file = atac_signal_file,
    bed4_file = bed4_file,
    signal_names = atac_names
  )

  # Mean ATAC across replicates
  ATAC[, ATAC_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = atac_names]

  # Merge signals into D
  D <- merge(
    D,
    ATAC[, .(id, ATAC_mean)],
    by = "id",
    all.x = TRUE
  )

  D <- merge(
    D,
    CHIP[, c("id", chip_names), with = FALSE],
    by = "id",
    all.x = TRUE
  )

  # -------------------------------------------------------
  # Save full merged table
  # -------------------------------------------------------

  fwrite(
    D,
    file.path(outDir, sprintf("review_%s_DMR_fUM_region_ChIP_ATAC_signal.tsv", dmrset)),
    sep = "\t"
  )

  # -------------------------------------------------------
  # Summary: ATAC median by region and fUM bin
  # -------------------------------------------------------

  ATAC_sum <- D[
    !is.na(fUM_bin) & !is.na(dmr_loc) & !is.na(ATAC_mean),
    .(
      nDMR = uniqueN(id),
      ATAC_median = median(ATAC_mean, na.rm = TRUE),
      ATAC_mean = mean(ATAC_mean, na.rm = TRUE)
    ),
    by = .(dmr_loc, fUM_bin)
  ][order(dmr_loc, fUM_bin)]

  fwrite(
    ATAC_sum,
    file.path(outDir, sprintf("review_%s_ATAC_signal_by_fUMbin_region.tsv", dmrset)),
    sep = "\t"
  )

  # -------------------------------------------------------
  # Summary: ChIP median by region, fUM bin and mark
  # -------------------------------------------------------

  L_chip <- melt(
    D,
    id.vars = c("id", "chr", "start", "end", "dmr_loc", "fUM", "fUM_bin", "ATAC_mean"),
    measure.vars = chip_names,
    variable.name = "mark",
    value.name = "ChIP_signal"
  )

  L_chip[, mark := factor(
    mark,
    levels = c("H3K4me3", "H3K36me3", "H3K27ac", "H3K9ac", "Pol2", "H3K9me2", "H3K27me3")
  )]

  CHIP_sum <- L_chip[
    !is.na(fUM_bin) & !is.na(dmr_loc) & !is.na(ChIP_signal),
    .(
      nDMR = uniqueN(id),
      ChIP_median = median(ChIP_signal, na.rm = TRUE),
      ChIP_mean = mean(ChIP_signal, na.rm = TRUE)
    ),
    by = .(dmr_loc, fUM_bin, mark)
  ][order(dmr_loc, mark, fUM_bin)]

  fwrite(
    CHIP_sum,
    file.path(outDir, sprintf("review_%s_ChIP_signal_by_fUMbin_region_mark.tsv", dmrset)),
    sep = "\t"
  )

  # =======================================================
  # Plot 1:
  # ChIP continuous signal by fUM bin and region
  # =======================================================

  p_chip <- ggplot(
    CHIP_sum,
    aes(x = fUM_bin, y = ChIP_median, color = mark, group = mark)
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.8) +
    facet_grid(. ~ dmr_loc, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.spacing = grid::unit(0.35, "lines")
    ) +
    labs(
      title = paste0(dmrset, " ChIP-seq signal over DMRs by fUM bin"),
      x = "fUM bin",
      y = "Median ChIP-seq signal over DMRs",
      color = "ChIP mark"
    )

  ggsave(
    file.path(outPlot, sprintf("review_%s_ChIP_signal_by_fUMbin_region_mark.pdf", dmrset)),
    p_chip,
    width = 24,
    height = 16,
    units = "cm"
  )

  # =======================================================
  # Plot 2:
  # ATAC continuous signal by fUM bin and region
  # =======================================================

  p_atac <- ggplot(
    ATAC_sum,
    aes(x = fUM_bin, y = ATAC_median, color = dmr_loc, group = dmr_loc)
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.0) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(
      title = paste0(dmrset, " ATAC-seq signal over DMRs by fUM bin"),
      x = "fUM bin",
      y = "Median ATAC-seq signal over DMRs",
      color = "DMR region"
    )

  ggsave(
    file.path(outPlot, sprintf("review_%s_ATAC_signal_by_fUMbin_region.pdf", dmrset)),
    p_atac,
    width = 22,
    height = 10,
    units = "cm"
  )

  # =======================================================
  # Plot 3:
  # Combined ChIP + ATAC continuous signal
  #
  # Because ChIP and ATAC are different scales, this plot uses
  # a secondary axis. ATAC is rescaled to the ChIP axis range.
  # This is useful for visual comparison of trends, not absolute
  # magnitude comparison.
  # =======================================================

  # Use global scale within each DMR set to keep the ATAC axis
  # consistent across regions and ChIP marks.
  chip_min <- min(CHIP_sum$ChIP_median, na.rm = TRUE)
  chip_max <- max(CHIP_sum$ChIP_median, na.rm = TRUE)

  atac_min <- min(ATAC_sum$ATAC_median, na.rm = TRUE)
  atac_max <- max(ATAC_sum$ATAC_median, na.rm = TRUE)

  if (
    !is.finite(chip_min) || !is.finite(chip_max) ||
    !is.finite(atac_min) || !is.finite(atac_max) ||
    chip_min == chip_max || atac_min == atac_max
  ) {
    warning("Skipping combined ChIP + ATAC plot for ", dmrset, " because signal range is invalid.")
  } else {

    scale_fac <- (chip_max - chip_min) / (atac_max - atac_min)

    ATAC_sum[, ATAC_scaled := chip_min + (ATAC_median - atac_min) * scale_fac]

    p_combined <- ggplot() +

      # ChIP signal
      geom_line(
        data = CHIP_sum,
        aes(x = fUM_bin, y = ChIP_median, color = mark, group = mark),
        linewidth = 0.6
      ) +
      geom_point(
        data = CHIP_sum,
        aes(x = fUM_bin, y = ChIP_median, color = mark, group = mark),
        size = 0.8
      ) +

      # ATAC signal
      geom_line(
        data = ATAC_sum,
        aes(x = fUM_bin, y = ATAC_scaled, group = 1),
        linewidth = 0.8,
        linetype = "dashed",
        color = "black"
      ) +
      geom_point(
        data = ATAC_sum,
        aes(x = fUM_bin, y = ATAC_scaled),
        size = 1.2,
        shape = 17,
        color = "black"
      ) +

      facet_grid(. ~ dmr_loc) +

      scale_y_continuous(
        name = "Median ChIP-seq signal over DMRs",
        sec.axis = sec_axis(
          transform = ~ (. - chip_min) / scale_fac + atac_min,
          name = "Median ATAC-seq signal over DMRs"
        )
      ) +

      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing = grid::unit(0.35, "lines")
      ) +

      labs(
        title = paste0(dmrset, " ChIP-seq and ATAC-seq signal over DMRs by fUM bin"),
        subtitle = "Solid colored lines: ChIP-seq marks; black dashed line: ATAC-seq",
        x = "fUM bin",
        color = "ChIP mark"
      )

    ggsave(
      file.path(outPlot, sprintf("review_%s_ChIP_ATAC_signal_by_fUMbin_region.pdf", dmrset)),
      p_combined,
      width = 24,
      height = 16,
      units = "cm"
    )
  }

  # =======================================================
  # Plot 4:
  # Normalized signal plot
  #
  # This is often safer than double y-axis for reviews.
  # It shows whether each signal increases/decreases across fUM bins.
  # =======================================================

  CHIP_norm <- copy(CHIP_sum)
  CHIP_norm[, signal_type := as.character(mark)]
  CHIP_norm[, signal := ChIP_median]
  CHIP_norm <- CHIP_norm[, .(dmr_loc, fUM_bin, signal_type, signal)]

  ATAC_norm <- copy(ATAC_sum)
  ATAC_norm[, signal_type := "ATAC"]
  ATAC_norm[, signal := ATAC_median]
  ATAC_norm <- ATAC_norm[, .(dmr_loc, fUM_bin, signal_type, signal)]

  L_norm <- rbindlist(
    list(CHIP_norm, ATAC_norm),
    use.names = TRUE,
    fill = TRUE
  )

  L_norm[
    ,
    signal_scaled_0_1 := rescale_to_range(signal, new_min = 0, new_max = 1),
    by = .(dmr_loc, signal_type)
  ]

  L_norm[, signal_type := factor(
    signal_type,
    levels = c("ATAC", "H3K4me3", "H3K36me3", "H3K27ac", "H3K9ac", "Pol2", "H3K9me2", "H3K27me3")
  )]

  p_norm <- ggplot(
    L_norm[!is.na(signal_scaled_0_1)],
    aes(x = fUM_bin, y = signal_scaled_0_1, color = signal_type, group = signal_type)
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.8) +
    facet_grid(. ~ dmr_loc) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.spacing = grid::unit(0.35, "lines")
    ) +
    labs(
      title = paste0(dmrset, " normalized ChIP-seq and ATAC-seq signal over DMRs"),
      subtitle = "Each signal is scaled from 0 to 1 within each DMR region",
      x = "fUM bin",
      y = "Scaled median signal",
      color = "Signal"
    )

  ggsave(
    file.path(outPlot, sprintf("review_%s_ChIP_ATAC_signal_scaled_by_fUMbin_region.pdf", dmrset)),
    p_norm,
    width = 24,
    height = 16,
    units = "cm"
  )

  message("Finished ", dmrset)
}
