#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

options(scipen = 999)

# Run from the repository root.
old_file <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/results/01.01_correlations-per-window_2kb.bed"
new_file <- "/mnt/disk2/vibanez/11_extended-data-analysis/Fig03_extended-data/extended-data-fig03a_01.01_correlations-per-window_3kb.bed"
outDir <- "/mnt/disk2/vibanez/11_extended-data-analysis/Fig03_extended-data"

table_file <- file.path(outDir,"extended-data-fig03a_positive-negative-correlation-ratio.tsv")
plot_file <- file.path(outDir,  "extended-data-fig03a_positive-negative-correlation-ratio.pdf")

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

for (input_file in c(old_file, new_file)) {
  if (!file.exists(input_file)) {
    stop("Missing input file: ", input_file)
  }
}

required_columns <- c(
  "geneID",
  "DMR",
  "pos",
  "geneType",
  "dmr",
  "corrResult",
  "corrDirection"
)

read_correlations <- function(path) {
  dat <- fread(
    path,
    sep = "\t",
    data.table = TRUE,
    fill = TRUE,
    check.names = FALSE,
    na.strings = "NA",
    nThread = 10
  )

  missing_columns <- setdiff(required_columns, names(dat))
  if (length(missing_columns)) {
    stop(
      "Missing columns in ", path, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat[, ..required_columns]
}

old_correlations <- read_correlations(old_file)
new_correlations <- read_correlations(new_file)

old_core <- copy(old_correlations)
old_core[, pos := as.integer(pos) + 15L]

new_outer <- new_correlations[
  (pos >= 0L & pos <= 14L) | (pos >= 59L & pos <= 73L)
]
new_outer[, pos := as.integer(pos)]

all_correlations <- rbindlist(
  list(old_core, new_outer),
  use.names = TRUE,
  fill = FALSE
)

all_correlations[, `:=`(
  geneID = as.character(geneID),
  DMR = as.character(DMR),
  geneType = as.character(geneType),
  dmr = as.character(dmr),
  corrResult = as.character(corrResult),
  corrDirection = as.character(corrDirection)
)]

all_correlations <- unique(
  all_correlations,
  by = required_columns
)

# Unfiltered analysis: retain every qualifying DMR–gene association.
plot_data <- all_correlations %>%
  filter(
    pos >= 5L,
    pos <= 68L,
    geneType == "gene",
    dmr %in% c("C-DMR", "CG-DMR"),
    corrResult == "correlated",
    corrDirection %in% c("positive", "negative")
  ) %>%
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
    # Preserve the display rule used for the original unfiltered panel.
    PositiveNegativeRatio = if_else(pos >= 40L, 1, RawRatio)
  )

fwrite(
  plot_data,
  table_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

x_breaks <- c(4.5, 14.5, 24.5, 34.5, 40, 48.5, 58.5, 68.5)
x_labels <- c("-3 kb", "-2 kb", "-1 kb", "TSS", "TES", "+1 kb", "+2 kb", "+3 kb")

p <- ggplot(
  plot_data,
  aes(
    x = pos,
    y = PositiveNegativeRatio,
    colour = dmr,
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
    colour = "grey30"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    colour = "grey40"
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(
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
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5)) +
  coord_cartesian(ylim = c(0, 1.5)) +
  labs(
    x = NULL,
    y = "Positive/negative correlation ratio",
    colour = "DMR type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(
  plot_file,
  p,
  width = 18,
  height = 12,
  units = "cm"
)

message("Saved table: ", table_file)
message("Saved plot: ", plot_file)
