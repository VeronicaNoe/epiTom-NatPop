# =========================================================
# Plot ChIP-seq mark overlap and ATAC-seq signal by fUM bin
# ATAC is plotted on a secondary y-axis
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
#DMRSETS <- c("CG-DMR")

# =========================================================
# Helper functions
# =========================================================

read_ids <- function(f) {
  if (!file.exists(f)) return(character())
  x <- fread(f, header = FALSE)
  if (nrow(x) == 0) return(character())
  as.character(x[[4]])
}

fmt_n <- function(x) {
  fifelse(
    x >= 1e6, paste0(round(x / 1e6, 1), "M"),
    fifelse(
      x >= 1e3, paste0(round(x / 1e3, 1), "k"),
      as.character(x)
    )
  )
}

# =========================================================
# Main loop
# =========================================================

for (dmrset in DMRSETS) {
  
  message("Processing: ", dmrset)
  annot_file <- sprintf(
    "%s/ac_annotate-DMRs-M82/dmr_mark_annotations/%s/%s.dmr_chromatin_annot.tsv.gz",
    OUTROOT, dmrset, dmrset
  )
  
  freq_file <- sprintf(
    "%s/ad_UM-frequency/master/%s.methfreq.tsv.gz",
    OUTROOT, dmrset
  )
  
  atac_tab <- sprintf(
    "%s/aa_prepare-public-data/master/%s.ATAC.mean_by_DMR.tsv",
    OUTROOT, dmrset
  )
  
  bed4_file <- file.path(
    OUTROOT,
    "aa_prepare-public-data/master",
    paste0(dmrset, ".M82.bed4")
  )
  
  
  A <- fread(annot_file)
  
  F <- fread(freq_file)
  setnames(
    F,
    c("id", "chr", "pos", "end", "n_um", "n_m", "n_obs", "fUM", "fMeth")
  )
  
  F[, id := as.character(id)]
  
  AT <- fread(atac_tab)
  setnames(AT, c("chr", "start", "end", "ATAC_rep1", "ATAC_rep2"))
  
  AT[, start := as.integer(start)]
  AT[, end   := as.integer(end)]
  
  AT[, c("ATAC_rep1", "ATAC_rep2") := lapply(.SD, as.numeric),
     .SDcols = c("ATAC_rep1", "ATAC_rep2")]
  
  AT[, ATAC_mean := rowMeans(.SD, na.rm = TRUE),
     .SDcols = c("ATAC_rep1", "ATAC_rep2")]
  
  # BED4 used to recover the DMR id for ATAC table
  B <- fread(bed4_file, header = FALSE)
  setnames(B, c("chr", "start", "end", "id"))
  
  B[, start := as.integer(start)]
  B[, end   := as.integer(end)]
  B[, id    := as.character(id)]
  
  AT_with_id <- merge(
    AT,
    B,
    by = c("chr", "start", "end"),
    all.x = TRUE
  )
  
  
  A[, id := as.character(id)]
  
  D <- merge(
    A,
    F[, .(id, n_obs, fUM, fMeth)],
    by = "id",
    all.x = TRUE
  )
  
  D <- merge(
    D,
    AT_with_id[, .(id, ATAC_mean)],
    by = "id",
    all.x = TRUE
  )
  
  
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
  
  loc_map <- data.table(
    id = D$id,
    dmr_loc = NA_character_
  )
  
  loc_map[id %in% ids_prom, dmr_loc := "promoter"]
  loc_map[is.na(dmr_loc) & id %in% ids_geneTE,    dmr_loc := "gene-TE"]
  loc_map[is.na(dmr_loc) & id %in% ids_gene_noTE, dmr_loc := "gene"]
  loc_map[is.na(dmr_loc) & id %in% ids_te,        dmr_loc := "TE"]
  loc_map[is.na(dmr_loc) & id %in% ids_inter,     dmr_loc := "intergenic"]
  
  D <- merge(D, loc_map, by = "id", all.x = TRUE)
  
  D[, dmr_loc := factor(
    dmr_loc,
    levels = c("promoter", "gene-TE", "gene", "TE", "intergenic")
  )]
  
  print(D[, .N, by = dmr_loc][order(-N)])
  
  
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
  
  
  Ntab_loc <- D[
    !is.na(fUM_bin),
    .(nDMR = uniqueN(id)),
    by = .(dmr_loc, fUM_bin)
  ][order(dmr_loc, fUM_bin)]
  
  fwrite(
    Ntab_loc,
    file.path(outDir, sprintf("review_%s_N_by_fUMbin_region.tsv", dmrset)),
    sep = "\t"
  )
  
  print(Ntab_loc)
  
  
  p_atac <- ggplot(
    D[!is.na(fUM_bin) & !is.na(ATAC_mean) & !is.na(dmr_loc)],
    aes(x = fUM_bin, y = ATAC_mean, color = dmr_loc, group = dmr_loc)
  ) +
    stat_summary(fun = median, geom = "line", linewidth = 0.6) +
    stat_summary(fun = median, geom = "point", size = 0.8) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(
      title = paste0(dmrset, " ATAC median signal by fUM bin"),
      x = "fUM bin (0/0 frequency)",
      y = "ATAC median signal",
      color = "Region"
    )
  
  ggsave(
    file.path(outPlot, sprintf("review_%s_ATAC_median_by_fUMbin.pdf", dmrset)),
    p_atac,
    width = 22,
    height = 10,
    units = "cm"
  )
  
  
  marks_any <- c(
    "H3K4me3_any", "H3K27ac_any", "H3K9ac_any",
    "H3K36me3_any", "H3K9me2_any", "Pol2_any",
    "H3K27me3_any"
  )
  
  marks_frac <- c(
    "H3K4me3_frac", "H3K27ac_frac", "H3K9ac_frac",
    "H3K36me3_frac", "H3K9me2_frac", "Pol2_frac",
    "H3K27me3_frac"
  )
  
  Lf <- melt(
    D,
    id.vars = c("id", "dmr_loc", "fUM", "fUM_bin", "n_obs"),
    measure.vars = intersect(marks_frac, names(D)),
    variable.name = "mark",
    value.name = "mark_frac"
  )
  
  Lf[, mark := sub("_frac$", "", mark)]
  
  S2 <- Lf[!is.na(fUM_bin) & !is.na(mark_frac)]
  
  S2[, mark_any := as.integer(mark_frac > 0)]
  
  tab_any <- S2[
    ,
    .(
      n = uniqueN(id),
      pct_any = 100 * mean(mark_any == 1)
    ),
    by = .(dmr_loc, mark, fUM_bin)
  ][order(dmr_loc, mark, fUM_bin)]
  
  
  p_any <- ggplot(
    tab_any,
    aes(x = fUM_bin, y = pct_any, color = dmr_loc, group = dmr_loc)
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.8) +
    facet_grid(. ~ mark) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(
      title = paste0(dmrset, " %DMRs with mark peak overlap by fUM bin"),
      x = "fUM bin (0/0 frequency)",
      y = "% DMRs with peak overlap",
      color = "Region"
    )
  
  ggsave(
    file.path(outPlot, sprintf("review_%s_pctAny_by_fUMbin_region.pdf", dmrset)),
    p_any,
    width = 28,
    height = 10,
    units = "cm"
  )
  
  
  mark_levels <- c(
    "H3K4me3", "H3K36me3", "H3K27ac", "H3K9ac",
    "Pol2", "H3K9me2", "H3K27me3"
  )
  
  tab_any[, mark := factor(
    mark,
    levels = intersect(mark_levels, unique(mark))
  )]
  
  p_any_loc <- ggplot(
    tab_any,
    aes(x = fUM_bin, y = pct_any, color = mark, group = mark)
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.8) +
    facet_grid(. ~ dmr_loc, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(
      title = paste0(dmrset, " %DMRs with peak overlap by fUM bin, per region"),
      x = "fUM bin (0/0 frequency)",
      y = "% DMRs with peak overlap",
      color = "Mark"
    )
  
  ggsave(
    file.path(outPlot, sprintf("review_%s_pctAny_by_fUMbin_byRegion_allMarks.pdf", dmrset)),
    p_any_loc,
    width = 24,
    height = 16,
    units = "cm"
  )
  
  
  ATAC_sum <- D[
    !is.na(fUM_bin) & !is.na(ATAC_mean) & !is.na(dmr_loc),
    .(
      ATAC_median = median(ATAC_mean, na.rm = TRUE),
      n_ATAC = uniqueN(id)
    ),
    by = .(dmr_loc, fUM_bin)
  ][order(dmr_loc, fUM_bin)]
  
  fwrite(
    ATAC_sum,
    file.path(outDir, sprintf("review_%s_ATAC_median_by_fUMbin_region.tsv", dmrset)),
    sep = "\t"
  )
  
  
  left_min <- 0
  left_max <- 100
  
  atac_min <- min(ATAC_sum$ATAC_median, na.rm = TRUE)
  atac_max <- max(ATAC_sum$ATAC_median, na.rm = TRUE)
  
  if (!is.finite(atac_min) || !is.finite(atac_max) || atac_min == atac_max) {
    stop("ATAC_median has invalid or constant values; cannot build secondary axis.")
  }
  
  scale_fac <- (left_max - left_min) / (atac_max - atac_min)
  
  ATAC_sum[, ATAC_scaled := left_min + (ATAC_median - atac_min) * scale_fac]
  
  p_any_loc_atac <- ggplot() +
    
    # ChIP-seq marks
    geom_line(
      data = tab_any,
      aes(x = fUM_bin, y = pct_any, color = mark, group = mark),
      linewidth = 0.6
    ) +
    geom_point(
      data = tab_any,
      aes(x = fUM_bin, y = pct_any, color = mark, group = mark),
      size = 0.8
    ) +
    
    # ATAC-seq median signal
    geom_line(
      data = ATAC_sum,
      aes(x = fUM_bin, y = ATAC_scaled, group = 1),
      linewidth = 0.7,
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
      name = "% DMRs with ChIP peak overlap",
      limits = c(0, 100),
      sec.axis = sec_axis(
        trans = ~ (. - left_min) / scale_fac + atac_min,
        name = "ATAC median signal"
      )
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.spacing = grid::unit(0.4, "lines")
    ) +
    
    labs(
      title = paste0(dmrset, " ChIP marks and ATAC signal by fUM bin"),
      x = "fUM bin (0/0 frequency)",
      color = "ChIP mark"
    )
  
  ggsave(
    file.path(outPlot, sprintf("review_%s_pctAny_ATAC_by_fUMbin_byRegion_allMarks.pdf", dmrset)),
    p_any_loc_atac,
    width = 24,
    height = 16,
    units = "cm"
  )
  
  
  N_mark <- S2[
    !is.na(fUM_bin),
    .(
      nDMR = uniqueN(id),
      n_marked = uniqueN(id[mark_any == 1]),
      pct_marked = 100 * uniqueN(id[mark_any == 1]) / uniqueN(id)
    ),
    by = .(dmr_loc, mark, fUM_bin)
  ][order(dmr_loc, mark, fUM_bin)]
  
  N_mark[, n_unmarked := nDMR - n_marked]
  N_mark[, pct_marked := 100 * n_marked / pmax(nDMR, 1)]
  N_mark[, pct_unmarked := 100 - pct_marked]
  
  fwrite(
    N_mark,
    file.path(outDir, sprintf("review_%s_N_marked_by_fUMbin_region_mark.tsv", dmrset)),
    sep = "\t"
  )
  
  
  Lcount <- melt(
    N_mark,
    id.vars = c("dmr_loc", "mark", "fUM_bin", "nDMR"),
    measure.vars = c("n_marked", "n_unmarked"),
    variable.name = "status",
    value.name = "n"
  )
  
  Lcount[, status := fifelse(status == "n_marked", "marked", "unmarked")]
  Lcount[, status := factor(status, levels = c("marked", "unmarked"))]
  
  Lcount[, pct := 100 * n / sum(n), by = .(dmr_loc, mark, fUM_bin)]
  
  THRESH_PCT <- 10
  
  Lcount[, label := fifelse(pct >= THRESH_PCT, fmt_n(n), "")]
  
  p100 <- ggplot(Lcount, aes(x = fUM_bin, y = n, fill = status)) +
    geom_col(position = "fill", width = 0.9) +
    geom_text(
      aes(label = label),
      position = position_fill(vjust = 0.5),
      size = 2.2,
      check_overlap = TRUE
    ) +
    facet_grid(dmr_loc ~ mark) +
    scale_y_continuous(labels = function(x) 100 * x) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.spacing = grid::unit(0.25, "lines")
    ) +
    labs(
      title = paste0(dmrset, " marked vs unmarked by fUM bin"),
      x = "fUM bin",
      y = "% of DMRs",
      fill = "Status"
    )
  
  ggsave(
    file.path(outPlot, sprintf("review_%s_marked_vs_unmarked_100pct_counts_by_fUMbin_region_mark.pdf", dmrset)),
    p100,
    width = 32,
    height = 16,
    units = "cm"
  )
  
  
  bed_out_dir <- file.path(OUTROOT, "deeptools_beds", dmrset)
  dir.create(bed_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (reg in levels(D$dmr_loc)) {
    for (b in levels(D$fUM_bin)) {
      
      sub <- D[
        dmr_loc == reg &
          fUM_bin == b &
          !is.na(chr) &
          !is.na(start) &
          !is.na(end),
        .(chr, start, end)
      ]
      
      b_str <- gsub("\\.", "p", b)
      
      out_bed <- file.path(
        bed_out_dir,
        sprintf("%s.%s.fUM_%s.bed", dmrset, reg, b_str)
      )
      
      fwrite(
        sub,
        out_bed,
        sep = "\t",
        col.names = FALSE
      )
    }
  }
  
  message("Finished: ", dmrset)
}
