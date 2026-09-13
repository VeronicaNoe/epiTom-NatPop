#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

options(scipen = 999)

# Run from the repository root.
input_file <- "10_data-analysis/Fig4/results/01.01_correlations-per-window_2kb.bed"
outDir <- "10_data-analysis/Fig4/results"
outPlot <- "10_data-analysis/Fig4/plots"

table_file <- file.path(outDir,"02.1_positive-negative-correlation-ratio_2kb.tsv")
plot_file <- file.path(outPlot,"02.1_positive-negative-correlation-ratio_2kb.pdf")

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
dir.create(outPlot, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
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

correlations <- fread(
  input_file,
  sep = "\t",
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = "NA",
  nThread = 10
)

missing_columns <- setdiff(required_columns, names(correlations))
if (length(missing_columns)) {
  stop(
    "Missing columns in ", input_file, ": ",
    paste(missing_columns, collapse = ", ")
  )
}

correlations <- correlations[, ..required_columns]
correlations[, `:=`(
  pos = as.integer(pos),
  geneID = as.character(geneID),
  DMR = as.character(DMR),
  geneType = as.character(geneType),
  dmr = as.character(dmr),
  corrResult = as.character(corrResult),
  corrDirection = as.character(corrDirection)
)]

correlations <- unique(correlations, by = required_columns)

# Native 2-kb coordinate system:
#   positions 0:19  = -2 kb to TSS
#   positions 20:23 = gene body
#   positions 24:43 = TES to +2 kb
plot_data <- correlations %>%
  filter(
    pos >= 0L,
    pos <= 43L,
    geneType == "gene",
    dmr %in% c("C-DMR", "CG-DMR"),
    corrResult == "correlated",
    corrDirection %in% c("positive", "negative")
  ) %>%
  distinct(pos, dmr, corrDirection, DMR) %>%
  count(pos, dmr, corrDirection, name = "count") %>%
  complete(
    pos = 0L:43L,
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
    # Equivalent to pos >= 40 after shifting the 2-kb coordinates by +15.
    PositiveNegativeRatio = if_else(pos >= 25L, 1, RawRatio)
  )

fwrite(
  plot_data,
  table_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

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
    xmin = 19.5,
    xmax = 25,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey85",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = c(19.5, 25),
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
    breaks = c(-0.5, 9.5, 19.5, 25, 33.5, 43.5),
    labels = c("-2 kb", "-1 kb", "TSS", "TES", "+1 kb", "+2 kb"),
    limits = c(-0.5, 43.5),
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
