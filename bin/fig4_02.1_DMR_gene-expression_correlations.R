suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

options(scipen = 999)

################################################################################
# Paths
################################################################################

in_file <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/results/01.01_correlations-per-window_3kb.bed"
all_genes_file <- "/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/allGenes.bed"
outPlot <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/"
outDir  <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/"

################################################################################
# Load the merged correlation table
#   pos  0:14 = -3.5 kb to -2 kb
#   pos 15:34 = -2 kb to TSS
#   pos 35:38 = gene body
#   pos 39:58 = TES to +2 kb
#   pos 59:73 = +2 kb to +3.5 kb
################################################################################

allCorrelations <- fread(
  in_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = "NA",
  nThread = 10
)

################################################################################
# Check required columns
################################################################################

required_columns <- c(
  "geneID",
  "DMR",
  "pos",
  "geneType",
  "dmr",
  "corrResult",
  "corrDirection"
)

missing_columns <- setdiff(required_columns, names(allCorrelations))

if (length(missing_columns) > 0L) {
  stop(
    "Missing required columns in ", in_file, ": ",
    paste(missing_columns, collapse = ", ")
  )
}

################################################################################
# Standardize types
################################################################################

allCorrelations[, `:=`(
  pos = as.integer(pos),
  geneID = as.character(geneID),
  DMR = as.character(DMR),
  geneType = as.character(geneType),
  dmr = as.character(dmr),
  corrResult = as.character(corrResult),
  corrDirection = as.character(corrDirection)
)]


allCorrelations <- unique(allCorrelations,by = required_columns)

################################################################################
# Load gene coordinates
################################################################################


allGenes <- fread(all_genes_file,sep = "\t",header = FALSE,col.names = c("chr","gene_start","gene_end","strand","geneID"),data.table = TRUE)

allGenes[,  `:=`(
    gene_start = as.integer(gene_start),
    gene_end = as.integer(gene_end),
    strand = as.character(strand),
    geneID = as.character(geneID))]

allGenes <- unique(allGenes,by = c(
    "chr",
    "gene_start",
    "gene_end",
    "strand",
    "geneID"  ))

# Use only focal genes represented in the correlation table
focalGenes <- allGenes[geneID %in% unique(allCorrelations$geneID)]


# All gene bodies used as potential neighboring genes

geneBodies <- allGenes[,.(chr,neighbor_start = gene_start,neighbor_end = gene_end - 1L,neighbor_gene = geneID)]

setkey(geneBodies,chr,neighbor_start,neighbor_end)


# Create exact 100-bp windows between 1 and 3.5 kb
#
# Upstream:
#   pos 0  = 3.4–3.5 kb upstream
#   pos 24 = 1.0–1.1 kb upstream
#
# Downstream:
#   pos 49 = 1.0–1.1 kb downstream
#   pos 73 = 3.4–3.5 kb downstream

positions_to_test <- data.table(pos = c(5L:34L, 39L:68L))

gene_position_windows <- focalGenes[,.(pos = positions_to_test$pos),by = .(chr,gene_start,gene_end,strand,geneID)]

gene_position_windows[,side := fifelse(pos <= 34L,"upstream","downstream")]

# Relative distance represented by each position
gene_position_windows[pos <= 34L,`:=`(distance_near = 3400L - pos * 100L,distance_far  = 3500L - pos * 100L)]

gene_position_windows[pos >= 39L,`:=`(distance_near = (pos - 39L) * 100L,distance_far  = (pos - 38L) * 100L  )]

# Convert relative positions to strand-aware genomic intervals


gene_position_windows[,window_start := fcase(side == "upstream" &strand == "+",gene_start - distance_far,
                                             side == "upstream" &strand == "-",gene_end + distance_near,
                                             side == "downstream" &strand == "+",gene_end + distance_near,
                                             side == "downstream" &strand == "-",gene_start - distance_far)]

gene_position_windows[,
  window_end := fcase(side == "upstream" & strand == "+", gene_start - distance_near - 1L,
                      side == "upstream" &strand == "-",gene_end + distance_far - 1L,
                      side == "downstream" &strand == "+",gene_end + distance_far - 1L,
                      side == "downstream" &strand == "-",gene_start - distance_near - 1L  )]

gene_position_windows[  window_start < 0L,  window_start := 0L]

gene_position_windows <- gene_position_windows[  window_end >= window_start]

gene_position_windows <- gene_position_windows[  ,  .(chr,focal_gene = geneID,pos,side,distance_near,distance_far,window_start,window_end)]


# Find another gene body overlapping each focal gene-position window


position_neighbor_hits <- foverlaps(x = gene_position_windows,y = geneBodies,by.x = c("chr","window_start","window_end"),by.y = c("chr","neighbor_start","neighbor_end"),type = "any",nomatch = 0L)

# Do not classify the focal gene itself as a neighboring gene
position_neighbor_hits <- position_neighbor_hits[focal_gene != neighbor_gene]

position_neighbor_hits[,overlap_bp := pmax(0L,pmin(window_end,neighbor_end) -pmax(window_start,neighbor_start) +1L)]

position_neighbor_hits <- unique(position_neighbor_hits,by = c("focal_gene","pos","neighbor_gene"))


# Collapse to one row per focal gene and position


position_neighbor_summary <- position_neighbor_hits[,.(other_gene_overlaps_position = TRUE,
                                                       n_other_genes_at_position =uniqueN(neighbor_gene),
                                                       other_genes_at_position =paste(sort(unique(neighbor_gene)),collapse = ","),
                                                       total_neighbor_overlap_bp =sum(overlap_bp)),
                                                    by = .(geneID = focal_gene,pos)]

fwrite(position_neighbor_summary,paste0(outDir,"01.00_gene-neighbors_by-position_1to3.5kb.tsv"),sep = "\t",quote = FALSE,na = "NA")


# Add positional overlap information to allCorrelations


allCorrelations[,`:=`(other_gene_overlaps_position = FALSE,
                      n_other_genes_at_position = 0L,
                      other_genes_at_position = NA_character_,
                      total_neighbor_overlap_bp = 0L)]

allCorrelations[position_neighbor_summary,on = .(geneID,pos),`:=`(
  other_gene_overlaps_position =i.other_gene_overlaps_position,
  n_other_genes_at_position = i.n_other_genes_at_position,
  other_genes_at_position = i.other_genes_at_position,
  total_neighbor_overlap_bp = i.total_neighbor_overlap_bp)]


# Save the complete annotated association table
fwrite(allCorrelations,paste0(outDir,"01.00_correlations_with-neighbor-gene_1to3kb.tsv"),sep = "\t",quote = FALSE,na = "NA")

# Summarize the effect of the neighboring-gene filter

neighbor_filter_summary <- allCorrelations %>%
  filter(pos >= 5L,pos <= 68L,geneType == "gene",dmr %in% c("C-DMR","CG-DMR")) %>%
  group_by(pos,dmr) %>%
  summarise(nDMRs_before =n_distinct(DMR),nDMRs_removed =n_distinct(DMR[other_gene_overlaps_position]),nDMRs_retained =n_distinct(DMR[!other_gene_overlaps_position]),percentage_removed =100 *nDMRs_removed /nDMRs_before,.groups = "drop")

fwrite(neighbor_filter_summary,paste0(outDir,"01.00_neighbor-filter_summary.tsv"),sep = "\t",quote = FALSE,na = "NA")



# Identify TE bodies overlapping each exact gene-position window

all_tes_file <- "/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/allTEs.bed"

rawTEs <- fread(all_tes_file,sep = "\t",header = FALSE,data.table = TRUE,fill = TRUE)

setnames(rawTEs,old = names(rawTEs)[1:3],new = c("chr","te_start","te_end"))

# Use column 4 as TE identifier when available
if (ncol(rawTEs) >= 4L) {
  setnames(rawTEs,old = names(rawTEs)[4],new = "teID")
} else {
  rawTEs[,teID := paste0(chr,":",te_start,"-",te_end)]
}

teBodies <- rawTEs[,.(chr = as.character(chr),
    te_start = as.integer(te_start),
    # BED end is normally exclusive;
    # convert it to an inclusive coordinate for foverlaps()
    te_end = as.integer(te_end) - 1L,
    teID = as.character(teID))]

teBodies <- teBodies[!is.na(chr) & !is.na(te_start) & !is.na(te_end) & te_end >= te_start]

teBodies <- unique(teBodies,by = c("chr","te_start","te_end","teID"))

setkey(teBodies,chr,te_start,te_end)

# Find TEs overlapping each exact 100-bp position

position_te_hits <- foverlaps(x = gene_position_windows,  y = teBodies,
                              by.x = c("chr","window_start","window_end"),
                              by.y = c("chr","te_start","te_end"),type = "any",nomatch = 0L)

position_te_hits[,overlap_bp := pmax(0L,pmin(window_end, te_end) -pmax(window_start, te_start) +1L)]

position_te_hits <- unique(position_te_hits,by = c("focal_gene","pos","teID"))

# Collapse to one row per focal gene and position

position_te_summary <- position_te_hits[,.(te_overlaps_position = TRUE,
    n_TEs_at_position =uniqueN(teID),
    TEs_at_position =paste(sort(unique(teID)),collapse = ","),total_TE_overlap_bp =sum(overlap_bp)),by = .(geneID = focal_gene,pos)]

fwrite(position_te_summary,paste0(outDir,"01.00_TE-overlap_by-gene-position_1to3.5kb.tsv"),sep = "\t",quote = FALSE,na = "NA")

# Add positional TE overlap information to allCorrelations

allCorrelations[,`:=`(te_overlaps_position = FALSE,n_TEs_at_position = 0L,TEs_at_position = NA_character_,total_TE_overlap_bp = 0L)]

allCorrelations[position_te_summary,
                on = .(geneID,pos),
                `:=`(te_overlaps_position =i.te_overlaps_position,
                     n_TEs_at_position =i.n_TEs_at_position,
                     TEs_at_position =i.TEs_at_position,
                     total_TE_overlap_bp =i.total_TE_overlap_bp)]

# Save complete position-annotated association table

fwrite(allCorrelations,paste0(outDir, "01.00_correlations_with-gene-TE-position-overlaps.tsv"),sep = "\t",quote = FALSE,na = "NA")


base_data <- allCorrelations %>%
  mutate(
    gene_overlap = coalesce(other_gene_overlaps_position, FALSE),
    TE_overlap = coalesce(te_overlaps_position, FALSE),
    position_overlap = gene_overlap | TE_overlap
  ) %>%
  filter(
    pos >= 5L,
    pos <= 68L,
    geneType == "gene",
    dmr %in% c("C-DMR", "CG-DMR"),
    corrResult == "correlated",
    corrDirection %in% c("positive", "negative")
  )

################################################################################
# Positive/negative ratio function
################################################################################

make_ratio <- function(data,analysis_name,remove_gene = FALSE,remove_TE = FALSE) {

  filtered_data <- data

  if (remove_gene) {
    filtered_data <- filtered_data %>%
      filter(!gene_overlap)
  }

  if (remove_TE) {
    filtered_data <- filtered_data %>%
      filter(!TE_overlap)
  }

  n_genes <- n_distinct(filtered_data$geneID)
  n_dmrs <- n_distinct(filtered_data$DMR)

  filtered_data %>%
    distinct(pos, dmr, corrDirection, DMR) %>%
    count(pos, dmr, corrDirection, name = "count") %>%
    complete(
      pos = 5L:68L,
      dmr = c("C-DMR", "CG-DMR"),
      corrDirection = c("positive", "negative"),
      fill = list(count = 0L)
    ) %>%
    pivot_wider(
      id_cols = c(pos, dmr),
      names_from = corrDirection,
      values_from = count,
      values_fill = list(count = 0L)
    ) %>%
    rename(
      Positive = positive,
      Negative = negative
    ) %>%
    mutate(
      support = Positive + Negative,
      Positive_plot = if_else(Positive == 0L, 1L, as.integer(Positive)),
      Negative_plot = if_else(Negative == 0L, 1L, as.integer(Negative)),
      RawRatio = Positive_plot / Negative_plot,

      # Keep the requested downstream display rule
      PositiveNegativeRatio = if_else(pos >= 40L, 1, RawRatio),

      analysis = analysis_name,
      nGenes = n_genes,
      nDMRs = n_dmrs
    )
}

# Define all position-specific filters

filter_specs <- tibble(
  analysis = c(
    "All associations",
    "No neighboring gene at position",
    "No TE at position",
    "No neighboring gene or TE at position"
  ),
  remove_gene = c(FALSE, TRUE, FALSE, TRUE),
  remove_TE = c(FALSE, FALSE, TRUE, TRUE)
)

# Calculate ratio profiles for all filters

ratio_list <- lapply(
  seq_len(nrow(filter_specs)),
  function(i) {
    make_ratio(
      data = base_data,
      analysis_name = filter_specs$analysis[i],
      remove_gene = filter_specs$remove_gene[i],
      remove_TE = filter_specs$remove_TE[i]
    )
  }
)

ratio_comparison <- bind_rows(ratio_list) %>%
  mutate(
    analysis = factor(
      analysis,
      levels = filter_specs$analysis
    )
  )

################################################################################
# Save combined and separate ratio tables
################################################################################

fwrite(ratio_comparison,  paste0(outDir, "01.00_ratio_position-filter-comparison.tsv"),  sep = "\t",  quote = FALSE,  na = "NA")

for (i in seq_len(nrow(filter_specs))) {

  analysis_name <- filter_specs$analysis[i]
  file_name <- gsub("[^A-Za-z0-9]+", "_", tolower(analysis_name))

  fwrite(
    ratio_comparison %>% filter(analysis == analysis_name),
    paste0(outDir, "01.00_ratio_", file_name, ".tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

# Plot settings

x_breaks <- c(
  4.5, 14.5, 24.5, 34.5,
  40, 48.5, 58.5, 68.5
)

x_labels <- c(
  "-3 kb", "-2 kb", "-1 kb", "TSS",
  "TES", "+1 kb", "+2 kb", "+3 kb"
)

make_ratio_plot <- function(plot_data) {

  ggplot(
    plot_data,
    aes(
      x = pos,
      y = PositiveNegativeRatio,
      color = dmr,
      group = dmr
    )
  ) +
    annotate(
      "rect",
      xmin = 34.5,
      xmax = 40,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey85",
      alpha = 0.5
    ) +
    geom_vline(
      xintercept = c(34.5, 40),
      linetype = "dashed",
      color = "grey30"
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dotted",
      color = "grey40"
    ) +
    geom_line(
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = c(
        "C-DMR" = "#820a86",
        "CG-DMR" = "#ffca7b"
      )
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(4.5, 68.5),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1, 1.5)
    ) +
    coord_cartesian(
      ylim = c(0, 1.5)
    ) +
    labs(
      x = NULL,
      y = "Positive/negative correlation ratio",
      color = "DMR type"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

# Save one ratio plot per filter

for (i in seq_len(nrow(filter_specs))) {

  analysis_name <- filter_specs$analysis[i]
  file_name <- gsub("[^A-Za-z0-9]+", "_", tolower(analysis_name))

  p_tmp <- make_ratio_plot(
    ratio_comparison %>%
      filter(analysis == analysis_name)
  )

  ggsave(paste0(outPlot, "01_ratio_", file_name, ".pdf"),    p_tmp,    width = 32,
    height = 16,    units = "cm"  )
}

################################################################################
# Four-panel ratio comparison
################################################################################

p_ratio_filters <- make_ratio_plot(ratio_comparison) +
  facet_wrap(
    ~ analysis,
    ncol = 2
  ) +
  theme(
    strip.text = element_text(face = "bold")
  )

ggsave(
  paste0(outPlot, "01_ratio_position-filter-comparison_4panels.pdf"),
  p_ratio_filters,
  width = 36,
  height = 26,
  units = "cm"
)

################################################################################
# Summarize the effect of each filter in the -3 to -2 kb interval
################################################################################

region_effect <- ratio_comparison %>%
  filter(pos >= 5L, pos <= 14L) %>%
  group_by(analysis, dmr) %>%
  summarise(
    Positive = sum(Positive),
    Negative = sum(Negative),
    support = Positive + Negative,
    ratio = Positive / Negative,
    .groups = "drop"
  )

region_baseline <- region_effect %>%
  filter(analysis == "All associations") %>%
  select(
    dmr,
    baseline_ratio = ratio
  )

region_effect <- region_effect %>%
  left_join(region_baseline, by = "dmr") %>%
  mutate(
    ratio_change_percent =
      100 * (ratio / baseline_ratio - 1)
  )

fwrite(
  region_effect,
  paste0(outDir, "01.00_ratio_filter-effects_3to2kb.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# Classify each unique DMR at each position

composition_DMR <- base_data %>%
  group_by(pos, dmr, corrDirection, DMR) %>%
  summarise(
    any_clean_association = any(!position_overlap),
    any_gene_overlap = any(gene_overlap),
    any_TE_overlap = any(TE_overlap),
    .groups = "drop"
  ) %>%
  mutate(
    position_assessed =
      pos %in% c(5L:34L, 39L:68L),

    overlap_class = case_when(
      !position_assessed ~
        "Gene body: retained",

      any_clean_association ~
        "Retained: no overlap",

      any_gene_overlap & any_TE_overlap ~
        "Removed: gene + TE",

      any_gene_overlap ~
        "Removed: gene",

      any_TE_overlap ~
        "Removed: TE",

      TRUE ~
        "Retained: no overlap"
    ),

    overlap_class = factor(
      overlap_class,
      levels = c(
        "Retained: no overlap",
        "Removed: TE",
        "Removed: gene",
        "Removed: gene + TE",
        "Gene body: retained"
      )
    )
  )

# Count and calculate percentages within positive and negative associations

composition_data <- composition_DMR %>%
  count(
    pos,
    dmr,
    corrDirection,
    overlap_class,
    name = "nDMR"
  ) %>%
  group_by(pos, dmr, corrDirection) %>%
  mutate(
    total_DMRs = sum(nDMR),
    percentage = 100 * nDMR / total_DMRs,
    signed_percentage = if_else(
      corrDirection == "negative",
      -percentage,
      percentage
    )
  ) %>%
  ungroup()

fwrite(
  composition_data,
  paste0(outDir, "01.00_position-overlap_composition_positive-negative.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


composition_difference <- composition_data %>%
  group_by(pos, dmr, overlap_class) %>%
  summarise(
    positive = sum(
      percentage[corrDirection == "positive"],
      na.rm = TRUE
    ),
    negative = sum(
      percentage[corrDirection == "negative"],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    negative_minus_positive =
      negative - positive
  )

fwrite(
  composition_difference,
  paste0(outDir, "01.00_position-overlap_negative-minus-positive.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Diverging stacked composition plot
#
# Positive associations are above zero.
# Negative associations are below zero.

p_composition <- ggplot(
  composition_data,
  aes(
    x = pos,
    y = signed_percentage,
    fill = overlap_class
  )
) +
  annotate(
    "rect",
    xmin = 34.5,
    xmax = 40,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey95",
    alpha = 0.5
  ) +
  geom_col(
    width = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  geom_vline(
    xintercept = c(34.5, 40),
    linetype = "dashed",
    color = "grey30"
  ) +
  scale_fill_manual(
    values = c(
      "Retained: no overlap" = "grey70",
      "Removed: TE" = "#E69F00",
      "Removed: gene" = "#0072B2",
      "Removed: gene + TE" = "#CC79A7",
      "Gene body: retained" = "grey92"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    limits = c(4.5, 68.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = c(-100, -50, 0, 50, 100),
    labels = function(x) paste0(abs(x), "%"),
    limits = c(-100, 100)
  ) +
  facet_grid(
    dmr ~ .
  ) +
  labs(
    x = NULL,
    y = "Share of correlated DMRs",
    fill = NULL,
    subtitle =
      "Positive correlations above zero; negative correlations below zero"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(
  paste0(outPlot, "02_position-overlap_composition_positive-negative.pdf"),
  p_composition,
  width = 36,
  height = 24,
  units = "cm"
)
