#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(mgcv)
})
# =============================================================================
# 0. Paths and parameters
# =============================================================================

setwd("10_data-analysis/Fig4/ad_snpEff")
full_genes_file <- "all_epialleles_plus_reference.genes.tsv"
full_variants_file <- "all_epialleles_plus_reference.variants.tsv"
enrichment_file <- "../results/03.02.accGroups_teM-frequency_PIM-SLC.tsv"
expression_file <- "../results/05.03.teM-geneExpression.tsv"
outDir <- "../results"
outPlot <- "../plots"
min_called <- 20L
fdr_cutoff <- 0.05
random_seed <- 1234L
set.seed(random_seed)
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
dir.create(outPlot, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(full_genes_file)) stop("Missing gene table: ", full_genes_file)
if (!file.exists(full_variants_file)) stop("Missing variant table: ", full_variants_file)
if (!file.exists(enrichment_file)) stop("Missing enrichment table: ", enrichment_file)
if (!file.exists(expression_file)) stop("Missing expression table: ", expression_file)
impact_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
strict_plof <- c("frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant", "start_lost", "transcript_ablation", "exon_loss_variant")
warning_terms <- c("WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS", "WARNING_TRANSCRIPT_NO_STOP_CODON", "WARNING_TRANSCRIPT_NO_START_CODON")
plof_class_levels <- c("No well-called pLoF", "Rare (<5%)", "Common (5-<50%)", "Frequent (50-<90%)", "Near-fixed (>=90%)")
class_colours <- c("Reference" = "#000000", "gbM" = "#ffca7b", "teM" = "#820a86", "PE" = "#4C78A8")
enrichment_colours <- c("Not enriched" = "grey70", "Enriched" = "#820a86")
expression_colours <- c("Downregulated" = "#d95f5f", "Upregulated" = "#4c78a8", "Non-DEG" = "grey65")
frequency_colours <- c("All genes with teM >0" = "black", "Domestication-enriched" = "#820a86", "Enriched downregulated" = "#d95f5f")

# =============================================================================
# 1. Helpers
# =============================================================================

normalize_gene <- function(x) x |> as.character() |> str_trim() |> str_remove("^(gene:|mRNA:)") |> str_remove("\\.[0-9]+$")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
bool <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}
repair_names <- function(x) {
  x[is.na(x) | x == ""] <- "unnamed"
  make.unique(x, sep = "_")
}
read_auto <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  reader <- if (str_count(first_line, "\t") >= str_count(first_line, ",")) read_tsv else read_csv
  reader(path, show_col_types = FALSE, progress = FALSE, name_repair = repair_names)
}
first_col <- function(dat, candidates, description, required = TRUE) {
  hit <- intersect(candidates, names(dat))
  if (length(hit)) return(hit[[1]])
  if (!required) return(NA_character_)
  stop("Cannot find ", description, ". Tried: ", paste(candidates, collapse = ", "), "\nAvailable columns: ", paste(names(dat), collapse = ", "))
}
format_p <- function(p) case_when(is.na(p) ~ "NA", p < 2.2e-16 ~ "<2.2e-16", p < 0.001 ~ formatC(p, format = "e", digits = 1), TRUE ~ formatC(p, format = "f", digits = 3))
n_distinct_if <- function(values, condition) n_distinct(values[replace_na(condition, FALSE)], na.rm = TRUE)
resolve_highest_impact <- function(x) map_chr(str_split(replace_na(as.character(x), ""), "[;&,]"), function(v) {
  v <- intersect(impact_levels, toupper(str_trim(v)))
  if (!length(v)) return(NA_character_)
  impact_levels[[min(match(v, impact_levels))]]
})
save_plot <- function(plot, filename, width, height) {
  ggsave(
    file.path(outPlot, filename),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}
binomial_summary <- function(dat, group_col, outcomes) {
  dat <- dat |> distinct(gene_id, .keep_all = TRUE)
  map_dfr(outcomes, function(outcome) dat |> filter(!is.na(.data[[group_col]]), !is.na(.data[[outcome]])) |> group_by(group = .data[[group_col]]) |> summarise(n_genes = n_distinct(gene_id), n_positive = sum(.data[[outcome]]), .groups = "drop") |> mutate(bt = map2(n_positive, n_genes, ~ binom.test(.x, .y)), percentage = 100 * n_positive / n_genes, CI95_low = 100 * map_dbl(bt, ~ .x$conf.int[[1]]), CI95_high = 100 * map_dbl(bt, ~ .x$conf.int[[2]]), outcome = outcome) |> select(-bt))
}
fisher_pair <- function(dat, group_col, exposed, reference, outcome, family) {
  x <- dat |> distinct(gene_id, .keep_all = TRUE) |> filter(.data[[group_col]] %in% c(exposed, reference), !is.na(.data[[outcome]]))
  n_exp <- sum(x[[group_col]] == exposed); n_ref <- sum(x[[group_col]] == reference)
  if (!n_exp || !n_ref) return(tibble(family = family, outcome = outcome, exposed_group = exposed, reference_group = reference, n_exposed = n_exp, positive_exposed = NA_integer_, pct_exposed = NA_real_, n_reference = n_ref, positive_reference = NA_integer_, pct_reference = NA_real_, odds_ratio = NA_real_, CI95_low = NA_real_, CI95_high = NA_real_, p_value = NA_real_))
  exp_pos <- sum(x[[group_col]] == exposed & x[[outcome]]); ref_pos <- sum(x[[group_col]] == reference & x[[outcome]])
  tab <- matrix(c(exp_pos, n_exp - exp_pos, ref_pos, n_ref - ref_pos), nrow = 2, byrow = TRUE)
  ft <- fisher.test(tab)
  tibble(family = family, outcome = outcome, exposed_group = exposed, reference_group = reference, n_exposed = n_exp, positive_exposed = exp_pos, pct_exposed = 100 * exp_pos / n_exp, n_reference = n_ref, positive_reference = ref_pos, pct_reference = 100 * ref_pos / n_ref, odds_ratio = unname(ft$estimate), CI95_low = ft$conf.int[[1]], CI95_high = ft$conf.int[[2]], p_value = ft$p.value)
}
pairwise_against_reference <- function(dat, group_col, groups, reference, outcomes, family) crossing(exposed_group = setdiff(groups, reference), outcome = outcomes) |> pmap_dfr(~ fisher_pair(dat, group_col, ..1, reference, ..2, family)) |> mutate(p_adjusted_BH = p.adjust(p_value, method = "BH"))
composition_test <- function(dat, group_col, status_col, family) {
  x <- dat |> distinct(gene_id, .keep_all = TRUE) |> filter(!is.na(.data[[group_col]]), !is.na(.data[[status_col]]))
  tab <- table(as.character(x[[group_col]]), as.character(x[[status_col]]))
  if (nrow(tab) < 2 || ncol(tab) < 2) return(tibble(family = family, method = NA_character_, statistic = NA_real_, degrees_freedom = NA_real_, p_value = NA_real_, min_expected = NA_real_, cramer_v = NA_real_, n_genes = sum(tab)))
  asym <- suppressWarnings(chisq.test(tab, correct = FALSE)); min_expected <- min(asym$expected)
  if (min_expected < 5) { tested <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 100000)); method <- "Pearson chi-square with Monte Carlo P value"; df <- NA_real_ } else { tested <- asym; method <- "Pearson chi-square"; df <- unname(asym$parameter) }
  n <- sum(tab); denom <- min(nrow(tab) - 1, ncol(tab) - 1)
  tibble(family = family, method = method, statistic = unname(asym$statistic), degrees_freedom = df, p_value = tested$p.value, min_expected = min_expected, cramer_v = ifelse(denom > 0, sqrt(unname(asym$statistic) / (n * denom)), NA_real_), n_genes = n)
}
class_summary <- function(dat, group_col) {
  dat |>
    mutate(group = .data[[group_col]]) |>
    filter(!is.na(group), !is.na(gene_pLoF_class)) |>
    distinct(gene_id, group, .keep_all = TRUE) |>
    count(group, pLoF_class = gene_pLoF_class, name = "n_genes") |>
    complete(group, pLoF_class = factor(plof_class_levels, levels = plof_class_levels), fill = list(n_genes = 0L)) |>
    group_by(group) |>
    mutate(total_genes = sum(n_genes), percentage = 100 * n_genes / total_genes) |>
    ungroup()
}
class_tests_against_reference <- function(dat, group_col, groups, reference, family) {
  crossing(
    exposed_group = setdiff(groups, reference),
    pLoF_class = plof_class_levels
  ) |>
    pmap_dfr(function(exposed_group, pLoF_class) {
      test_dat <- dat |>
        mutate(class_outcome = gene_pLoF_class == pLoF_class)
      fisher_pair(
        test_dat,
        group_col,
        exposed_group,
        reference,
        "class_outcome",
        family
      ) |>
        mutate(pLoF_class = pLoF_class)
    }) |>
    mutate(p_adjusted_BH = p.adjust(p_value, method = "BH"))
}

make_binary_count_plot <- function(summary_dat, tests, group_order, outcome_labels,
                                   colours, reference, title, subtitle, filename,
                                   width = 10, height = 7) {
  pdat <- summary_dat |>
    left_join(
      tests |>
        transmute(outcome, group = exposed_group, p_adjusted_BH),
      by = c("outcome", "group")
    ) |>
    mutate(
      group = factor(group, levels = group_order),
      outcome = factor(outcome, levels = names(outcome_labels), labels = outcome_labels),
      count_label = comma(n_positive),
      test_label = case_when(
        as.character(group) == reference ~ "ref.",
        is.na(p_adjusted_BH) ~ "",
        TRUE ~ paste0("q=", format_p(p_adjusted_BH))
      )
    )

  y_top <- max(pdat$n_positive, na.rm = TRUE)
  y_pad <- max(1, y_top * 0.035)
  y_limit <- max(1, y_top + 3.6 * y_pad)

  pdat <- pdat |>
    mutate(
      count_y = n_positive + y_pad,
      test_y = n_positive + 2.25 * y_pad
    )

  p <- ggplot(pdat, aes(group, n_positive, fill = group)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(aes(y = count_y, label = count_label), size = 2.7, vjust = 0) +
    geom_text(aes(y = test_y, label = test_label), size = 2.5, vjust = 0, fontface = "italic") +
    facet_wrap(~ outcome, scales = "fixed") +
    scale_fill_manual(values = colours, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, y_limit),
      labels = label_comma(),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Number of genes",
      title = title,
      subtitle = subtitle,
      caption = paste0(
        "Bars show gene counts. q values are BH-adjusted Fisher exact tests ",
        "of class prevalence versus the reference group; tests therefore account ",
        "for the different total number of genes in each group."
      )
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 28, hjust = 1),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 8, hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

  save_plot(p, filename, width, height)
  p
}

make_class_count_facets <- function(summary_dat, tests, group_order, group_labels,
                                    colours, reference, title, subtitle, filename,
                                    width = 11, height = 8, show_tests = TRUE) {
  pdat <- summary_dat |>
    mutate(
      group = as.character(group),
      pLoF_class = as.character(pLoF_class)
    )

  if (show_tests) {
    pdat <- pdat |>
      left_join(
        tests |>
          select(exposed_group, pLoF_class, p_adjusted_BH),
        by = c("group" = "exposed_group", "pLoF_class")
      )
  } else {
    pdat$p_adjusted_BH <- NA_real_
  }

  totals <- pdat |>
    distinct(group, total_genes) |>
    mutate(
      display_group = unname(group_labels[group]),
      axis_label = paste0(display_group, "\nN=", comma(total_genes))
    )
  axis_labels <- setNames(totals$axis_label, totals$group)

  pdat <- pdat |>
    mutate(
      group = factor(group, levels = group_order),
      pLoF_class = factor(pLoF_class, levels = plof_class_levels),
      count_label = comma(n_genes),
      test_label = case_when(
        !show_tests ~ "",
        as.character(group) == reference ~ "ref.",
        is.na(p_adjusted_BH) ~ "",
        TRUE ~ paste0("q=", format_p(p_adjusted_BH))
      )
    )

  y_top <- max(pdat$n_genes, na.rm = TRUE)
  y_pad <- max(1, y_top * 0.025)
  y_limit <- max(1, y_top + if (show_tests) 4.1 * y_pad else 2.5 * y_pad)

  pdat <- pdat |>
    mutate(
      count_y = n_genes + y_pad,
      test_y = n_genes + 2.25 * y_pad
    )

  p <- ggplot(pdat, aes(group, n_genes, fill = group)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(aes(y = count_y, label = count_label), size = 2.75, vjust = 0)

  if (show_tests) {
    p <- p +
      geom_text(
        aes(y = test_y, label = test_label),
        size = 2.5,
        vjust = 0,
        fontface = "italic"
      )
  }

  p <- p +
    facet_wrap(
      ~ pLoF_class,
      ncol = 3,
      scales = "fixed",
      labeller = as_labeller(plof_facet_labels)
    ) +
    scale_x_discrete(labels = axis_labels) +
    scale_fill_manual(values = colours, drop = FALSE, guide = "none") +
    scale_y_continuous(
      limits = c(0, y_limit),
      labels = label_comma(),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Number of genes",
      title = title,
      subtitle = subtitle,
      caption = if (show_tests) {
        paste0(
          "Numbers above bars are gene counts. q values are BH-adjusted Fisher ",
          "exact tests of class prevalence versus the reference group"
        )
      } else {
        paste0(
          "Numbers above bars are gene counts. Displayed subsets overlap, so no ",
          "independent group-comparison test is shown"
        )
      }
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      panel.spacing = unit(1.0, "lines"),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 8, hjust = 0),
      plot.margin = margin(8, 14, 8, 8)
    )

  save_plot(p, filename, width, height)
  p
}

model_terms <- function(model, model_name) {
  z <- summary(model)$coefficients; p_col <- grep("^Pr\\(", colnames(z), value = TRUE)[[1]]
  tibble(model = model_name, term = rownames(z), estimate = z[, "Estimate"], standard_error = z[, "Std. Error"], exponentiated_estimate = exp(z[, "Estimate"]), CI95_low = exp(z[, "Estimate"] - 1.96 * z[, "Std. Error"]), CI95_high = exp(z[, "Estimate"] + 1.96 * z[, "Std. Error"]), p_value = z[, p_col])
}
bootstrap_mean_ci <- function(x, B = 5000L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(tibble(mean = NA_real_, CI95_low = NA_real_, CI95_high = NA_real_))
  if (length(x) == 1) return(tibble(mean = x, CI95_low = x, CI95_high = x))
  b <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  tibble(mean = mean(x), CI95_low = unname(quantile(b, 0.025)), CI95_high = unname(quantile(b, 0.975)))
}

# =============================================================================
# 2. Read combined gene and variant tables
# =============================================================================

genes_raw <- read_auto(full_genes_file)
gene_col <- first_col(genes_raw, c("geneName", "gene", "gene_id", "GeneName"), "gene identifier")
genes <- genes_raw |>
  mutate(gene_id = normalize_gene(.data[[gene_col]])) |>
  filter(!is.na(gene_id), gene_id != "")
if ("n_unique_variants" %in% names(genes)) {
  genes <- genes |> arrange(gene_id, desc(num(n_unique_variants))) |> distinct(gene_id, .keep_all = TRUE)
} else {
  genes <- genes |> distinct(gene_id, .keep_all = TRUE)
}
required_gene_cols <- c("analysis_gene_class", "is_random_reference")
if (length(setdiff(required_gene_cols, names(genes)))) stop("Combined gene table is missing: ", paste(setdiff(required_gene_cols, names(genes)), collapse = ", "))
tem_frequency_col <- first_col(genes, c("epifreq_teM", "teM_frequency", "teM_freq"), "population teM frequency", required = FALSE)
genes <- genes |> mutate(is_random_reference = bool(is_random_reference), has_epiallele_data = if ("has_epiallele_data" %in% names(genes)) bool(has_epiallele_data) else !is_random_reference, has_snpeff_annotation = if ("has_snpeff_annotation" %in% names(genes)) bool(has_snpeff_annotation) else replace_na(num(n_unique_variants), 0) > 0, reference_warning_evaluable = if ("reference_warning_evaluable" %in% names(genes)) bool(reference_warning_evaluable) else has_snpeff_annotation, global_gene_class = if ("geneEpiallele" %in% names(genes)) as.character(geneEpiallele) else as.character(analysis_gene_class), comparison_gene_class = case_when(is_random_reference ~ "Reference", analysis_gene_class == "gbM" ~ "gbM", analysis_gene_class == "teM" ~ "teM", analysis_gene_class == "PE" ~ "PE", analysis_gene_class == "UM" ~ "UM", TRUE ~ "Other"), teM_frequency = if (!is.na(tem_frequency_col)) num(.data[[tem_frequency_col]]) else NA_real_)
if (any(!is.na(genes$teM_frequency)) && max(genes$teM_frequency, na.rm = TRUE) > 1.5) genes$teM_frequency <- genes$teM_frequency / 100
genes$annotated_pseudogene <- if ("annotated_pseudogene" %in% names(genes)) bool(genes$annotated_pseudogene) else if ("biotype" %in% names(genes)) str_detect(str_to_lower(replace_na(as.character(genes$biotype), "")), "pseudogene") else FALSE
warning_col <- first_col(genes, c("reference_warning_types", "reference_model_warning_types", "annotation_warnings", "transcript_warnings", "warnings", "warning"), "reference warning text", required = FALSE)
genes$reference_warning_types <- if (!is.na(warning_col)) as.character(genes[[warning_col]]) else NA_character_
warning_flag <- rep(FALSE, nrow(genes))
if ("reference_model_warning" %in% names(genes)) warning_flag <- warning_flag | bool(genes$reference_model_warning)
if ("problematic_transcript_model" %in% names(genes)) warning_flag <- warning_flag | bool(genes$problematic_transcript_model)
warning_flag <- warning_flag | str_detect(replace_na(genes$reference_warning_types, ""), str_c(warning_terms, collapse = "|"))
genes$reference_model_warning <- warning_flag
genes$reference_model_warning_test <- if_else(genes$reference_warning_evaluable, genes$reference_model_warning, NA)
variants_raw <- read_auto(full_variants_file)
v_gene <- first_col(variants_raw, c("geneName", "gene", "gene_id", "GeneName"), "gene identifier in variant table")
v_chr <- first_col(variants_raw, c("chromosome", "CHROM", "chr"), "chromosome")
v_pos <- first_col(variants_raw, c("position", "POS", "pos"), "position")
v_alt <- first_col(variants_raw, c("effect_allele", "ALT", "alt"), "effect allele")
v_effect <- first_col(variants_raw, c("effects", "effect", "plof_effects"), "effect field")
v_impact <- first_col(variants_raw, c("highest_impact", "impact", "impacts"), "impact field")
v_id <- first_col(variants_raw, c("variant_id", "ID"), "variant identifier", required = FALSE)
variants <- variants_raw |> mutate(gene_id = normalize_gene(.data[[v_gene]]), variant_id = if (!is.na(v_id)) as.character(.data[[v_id]]) else paste(.data[[v_chr]], .data[[v_pos]], .data[[v_alt]], sep = ":"), effects_value = as.character(.data[[v_effect]]), impact = resolve_highest_impact(.data[[v_impact]]), n_called = num(n_called), allele_frequency = num(allele_frequency), carried = !is.na(allele_frequency) & allele_frequency > 0, enough_calls = !is.na(n_called) & n_called >= min_called) |> filter(!is.na(gene_id), gene_id != "")
variants$is_pLoF <- if ("is_plof" %in% names(variants)) bool(variants$is_plof) else map_lgl(str_split(replace_na(variants$effects_value, ""), "[;&,]"), ~ any(str_trim(.x) %in% strict_plof))
variants$pLoF_effects_value <- if ("plof_effects" %in% names(variants)) as.character(variants$plof_effects) else map_chr(str_split(replace_na(variants$effects_value, ""), "[;&,]"), ~ paste(intersect(str_trim(.x), strict_plof), collapse = ";"))
variants <- variants |> arrange(gene_id, variant_id) |> distinct(gene_id, variant_id, .keep_all = TRUE)

# =============================================================================
# 3. One gene-level pLoF definition used everywhere
# =============================================================================

variant_counts <- variants |> group_by(gene_id) |> summarise(n_unique_variants = n_distinct(variant_id), n_carried_variants = n_distinct_if(variant_id, carried), n_carried_pLoF = n_distinct_if(variant_id, carried & is_pLoF), n_well_called_pLoF = n_distinct_if(variant_id, carried & enough_calls & is_pLoF), .groups = "drop")
gene_plof <- variants |> filter(is_pLoF, carried, enough_calls, !is.na(allele_frequency)) |> group_by(gene_id) |> summarise(gene_pLoF_AF = max(allele_frequency), median_pLoF_AF = median(allele_frequency), mean_pLoF_AF = mean(allele_frequency), .groups = "drop")
impact_gene <- variants |> filter(carried, enough_calls, !is.na(impact)) |> group_by(gene_id, impact) |> summarise(common_carried = any(allele_frequency >= 0.05), .groups = "drop") |> pivot_wider(names_from = impact, values_from = common_carried, values_fill = FALSE, names_glue = "common_carried_{impact}")
for (impact_name in impact_levels) if (!paste0("common_carried_", impact_name) %in% names(impact_gene)) impact_gene[[paste0("common_carried_", impact_name)]] <- FALSE
recalc <- c(names(variant_counts)[-1], "gene_pLoF_AF", "median_pLoF_AF", "mean_pLoF_AF", paste0("common_carried_", impact_levels), "gene_pLoF_class", "has_pLoF_rare", "has_pLoF_common", "has_pLoF_frequent", "has_pLoF_near_fixed", "any_well_called_pLoF")
gene_analysis <- genes |> select(-any_of(recalc)) |> left_join(variant_counts, by = "gene_id") |> left_join(gene_plof, by = "gene_id") |> left_join(impact_gene, by = "gene_id") |> mutate(across(any_of(c(names(variant_counts)[-1])), ~ as.integer(replace_na(num(.x), 0))), across(any_of(paste0("common_carried_", impact_levels)), ~ replace_na(as.logical(.x), FALSE)), gene_pLoF_AF = replace_na(gene_pLoF_AF, 0), median_pLoF_AF = replace_na(median_pLoF_AF, 0), mean_pLoF_AF = replace_na(mean_pLoF_AF, 0), gene_pLoF_class = case_when(gene_pLoF_AF == 0 ~ "No well-called pLoF", gene_pLoF_AF < 0.05 ~ "Rare (<5%)", gene_pLoF_AF < 0.50 ~ "Common (5-<50%)", gene_pLoF_AF < 0.90 ~ "Frequent (50-<90%)", TRUE ~ "Near-fixed (>=90%)"), gene_pLoF_class = factor(gene_pLoF_class, levels = plof_class_levels), has_pLoF_rare = gene_pLoF_class == "Rare (<5%)", has_pLoF_common = gene_pLoF_class == "Common (5-<50%)", has_pLoF_frequent = gene_pLoF_class == "Frequent (50-<90%)", has_pLoF_near_fixed = gene_pLoF_class == "Near-fixed (>=90%)", any_well_called_pLoF = gene_pLoF_AF > 0)
# =============================================================================
# 4. Reference, enrichment and expression data
# =============================================================================

reference_sample <- gene_analysis |> filter(is_random_reference) |> distinct(gene_id, .keep_all = TRUE)
if (!nrow(reference_sample)) stop("No random-reference genes found")
reference_n_actual <- n_distinct(reference_sample$gene_id)
global_groups <- c("Reference", "gbM", "teM", "PE")
global_dat <- bind_rows(reference_sample |> mutate(comparison_group = "Reference"), gene_analysis |> filter(comparison_gene_class %in% c("gbM", "teM", "PE")) |> mutate(comparison_group = comparison_gene_class)) |> distinct(gene_id, .keep_all = TRUE) |> mutate(comparison_group = factor(comparison_group, levels = global_groups))
enrich_raw <- read_auto(enrichment_file)
required_enrich <- c("geneName", "p.adjusted_fdr", "PIM_count", "SLC_count")
if (length(setdiff(required_enrich, names(enrich_raw)))) stop("Missing enrichment columns: ", paste(setdiff(required_enrich, names(enrich_raw)), collapse = ", "))
enrichment <- enrich_raw |> mutate(gene_id = normalize_gene(geneName), enrichment_padj = num(p.adjusted_fdr), enrichment_status = case_when(is.na(enrichment_padj) ~ NA_character_, enrichment_padj <= fdr_cutoff ~ "Enriched", TRUE ~ "Not enriched"), enrichment_direction = case_when(enrichment_status == "Enriched" & SLC_count > PIM_count ~ "SLC enriched", enrichment_status == "Enriched" & PIM_count > SLC_count ~ "PIM enriched", enrichment_status == "Enriched" ~ "Enriched, equal counts", enrichment_status == "Not enriched" ~ "Not enriched", TRUE ~ NA_character_)) |> arrange(gene_id, enrichment_padj) |> distinct(gene_id, .keep_all = TRUE) |> select(gene_id, enrichment_status, enrichment_direction, enrichment_padj, PIM_count, SLC_count)
expr_raw <- read_auto(expression_file)
expr_gene <- first_col(expr_raw, c("geneName", "gene", "gene_id", "GeneName"), "gene identifier in expression table")
expr_status <- first_col(expr_raw, c("significance", "significant_direction", "expression_group"), "expression status")
expr_p <- first_col(expr_raw, c("p.value", "p_value", "pvalue", "padj"), "expression P value", required = FALSE)
expr_logfc <- first_col(expr_raw, c("logFC", "log2FoldChange", "log2FC"), "expression log fold change", required = FALSE)
expression <- expr_raw |> mutate(gene_id = normalize_gene(.data[[expr_gene]]), expression_group_raw = as.character(.data[[expr_status]]), expression_group = case_when(str_to_lower(expression_group_raw) == "downregulated" ~ "Downregulated", str_to_lower(expression_group_raw) == "upregulated" ~ "Upregulated", str_to_lower(expression_group_raw) %in% c("neutral", "moderated", "non-deg", "non_deg") ~ "Non-DEG", TRUE ~ NA_character_), DE_status = case_when(expression_group %in% c("Downregulated", "Upregulated") ~ "DEG", expression_group == "Non-DEG" ~ "non-DEG", TRUE ~ NA_character_), expression_p_value = if (!is.na(expr_p)) num(.data[[expr_p]]) else NA_real_, expression_logFC = if (!is.na(expr_logfc)) num(.data[[expr_logfc]]) else NA_real_) |> arrange(gene_id, expression_p_value) |> distinct(gene_id, .keep_all = TRUE) |> select(gene_id, expression_group_raw, expression_group, DE_status, expression_p_value, expression_logFC)
tem_analysis <- enrichment |> left_join(gene_analysis |> mutate(in_gene_summary = TRUE), by = "gene_id") |> left_join(expression, by = "gene_id") |> mutate(in_gene_summary = replace_na(in_gene_summary, FALSE), expression_group = replace_na(expression_group, "Not classified"), DE_status = replace_na(DE_status, "Not classified"))
eligible_tem <- tem_analysis |> filter(in_gene_summary)
# =============================================================================
# Supplementary Figure 15a: pLoF-class counts by gene methylation class
# =============================================================================

reference_plot_label <- "Random reference"

global_group_labels <- c(
  "Reference" = reference_plot_label,
  "gbM" = "gbM",
  "teM" = "teM",
  "PE" = "PE"
)

plof_facet_labels <- c(
  "No well-called pLoF" = "No pLoF",
  "Rare (<5%)" = "Rare pLoF (<5%)",
  "Common (5-<50%)" = "Common pLoF (5–<50%)",
  "Frequent (50-<90%)" = "Frequent pLoF (50–<90%)",
  "Near-fixed (>=90%)" = "Near-fixed pLoF (≥90%)"
)

global_class_summary <- class_summary(global_dat, "comparison_group")
global_class_test <- composition_test(
  global_dat,
  "comparison_group",
  "gene_pLoF_class",
  "Global pLoF-class composition"
)
global_class_tests <- class_tests_against_reference(
  global_dat,
  "comparison_group",
  global_groups,
  "Reference",
  "Global pLoF classes versus random reference"
)

write_tsv(global_class_summary, file.path(outDir, "suppFig-15a_pLoF-class-composition.tsv"))
write_tsv(global_class_test, file.path(outDir, "suppFig-15a_pLoF-class-omnibus-test.tsv"))
write_tsv(global_class_tests, file.path(outDir, "suppFig-15a_pLoF-class-vs-reference-tests.tsv"))

p1 <- make_class_count_facets(
  summary_dat = global_class_summary,
  tests = global_class_tests,
  group_order = global_groups,
  group_labels = global_group_labels,
  colours = class_colours,
  reference = "Reference",
  title = "Strict pLoF frequency classes across gene methylation classes",
  subtitle = paste0(
    "Each gene is assigned once using the maximum AF of its well-called strict pLoF variants. ",
    global_class_test$method[[1]], ": P=", format_p(global_class_test$p_value[[1]]),
    "; Cramer's V=", sprintf("%.3f", global_class_test$cramer_v[[1]]), "."
  ),
  filename = "suppFig-15a_pLoF-class-counts.pdf",
  width = 11.5,
  height = 8.3,
  show_tests = TRUE
)

# Supplementary Figure 15c: pLoF-class counts by domestication enrichment
# =============================================================================
enrichment_groups <- c("Not enriched", "Enriched")
enrichment_group_labels <- c("Not enriched" = "Not enriched", "Enriched" = "Enriched")
enrichment_dat <- eligible_tem |>
  filter(enrichment_status %in% enrichment_groups) |>
  mutate(enrichment_status = factor(enrichment_status, levels = enrichment_groups))

enrichment_class_summary <- class_summary(enrichment_dat, "enrichment_status")
enrichment_class_test <- composition_test(
  enrichment_dat,
  "enrichment_status",
  "gene_pLoF_class",
  "Enrichment pLoF-class composition"
)
enrichment_class_tests <- class_tests_against_reference(
  enrichment_dat,
  "enrichment_status",
  enrichment_groups,
  "Not enriched",
  "Enriched versus non-enriched pLoF classes"
)
write_tsv(enrichment_class_summary, file.path(outDir, "suppFig-15c_pLoF-class-composition.tsv"))
write_tsv(enrichment_class_test, file.path(outDir, "suppFig-15c_pLoF-class-omnibus-test.tsv"))
write_tsv(enrichment_class_tests, file.path(outDir, "suppFig-15c_pLoF-class-vs-nonenriched-tests.tsv"))

make_class_count_facets(
  summary_dat = enrichment_class_summary,
  tests = enrichment_class_tests,
  group_order = enrichment_groups,
  group_labels = enrichment_group_labels,
  colours = enrichment_colours,
  reference = "Not enriched",
  title = "Strict pLoF-class counts in enriched and non-enriched teM genes",
  subtitle = paste0(
    enrichment_class_test$method[[1]], ": P=", format_p(enrichment_class_test$p_value[[1]]),
    "; Cramer's V=", sprintf("%.3f", enrichment_class_test$cramer_v[[1]]), "."
  ),
  filename = "suppFig-15c_pLoF-class-counts.pdf",
  width = 10.5,
  height = 7.5,
  show_tests = TRUE
)

# =============================================================================
# Supplementary Figure 15d: pLoF-class counts by expression group
# =============================================================================
expression_groups <- c("Downregulated", "Upregulated", "Non-DEG")
expression_group_labels <- c(
  "Downregulated" = "Downregulated",
  "Upregulated" = "Upregulated",
  "Non-DEG" = "Non-DEG"
)
expression_dat <- eligible_tem |>
  filter(enrichment_status == "Enriched", expression_group %in% expression_groups) |>
  mutate(expression_group = factor(expression_group, levels = expression_groups))

expression_class_summary <- class_summary(expression_dat, "expression_group")
expression_class_test <- composition_test(
  expression_dat,
  "expression_group",
  "gene_pLoF_class",
  "Expression pLoF-class composition"
)
expression_class_tests <- class_tests_against_reference(
  expression_dat,
  "expression_group",
  expression_groups,
  "Non-DEG",
  "Expression pLoF classes versus non-DEG"
)
write_tsv(expression_class_summary, file.path(outDir, "suppFig-15d_pLoF-class-composition.tsv"))
write_tsv(expression_class_test, file.path(outDir, "suppFig-15d_pLoF-class-omnibus-test.tsv"))
write_tsv(expression_class_tests, file.path(outDir, "suppFig-15d_pLoF-class-vs-non-DEG-tests.tsv"))

make_class_count_facets(
  summary_dat = expression_class_summary,
  tests = expression_class_tests,
  group_order = expression_groups,
  group_labels = expression_group_labels,
  colours = expression_colours,
  reference = "Non-DEG",
  title = "Strict pLoF-class counts among enriched teM expression groups",
  subtitle = paste0(
    expression_class_test$method[[1]], ": P=", format_p(expression_class_test$p_value[[1]]),
    "; Cramer's V=", sprintf("%.3f", expression_class_test$cramer_v[[1]]), "."
  ),
  filename = "suppFig-15d_pLoF-class-counts.pdf",
  width = 10.8,
  height = 7.5,
  show_tests = TRUE
)

# =============================================================================
# Supplementary Figure 15b: pLoF frequency versus teM frequency
# =============================================================================

continuous_dat <- gene_analysis |>
  left_join(enrichment |> select(gene_id, enrichment_status), by = "gene_id") |>
  left_join(expression |> select(gene_id, expression_group), by = "gene_id") |>
  mutate(
    teM_frequency_10pct = teM_frequency / 0.10,
    log_total_variants = log1p(n_unique_variants)
  ) |>
  distinct(gene_id, .keep_all = TRUE)

teM_continuous <- continuous_dat |>
  filter(
    !is_random_reference,
    is.finite(teM_frequency),
    teM_frequency > 0,
    teM_frequency <= 1,
    is.finite(gene_pLoF_AF)
  )

all_teM_positive <- teM_continuous |>
  mutate(analysis_group = "All genes with teM >0")

domestication_enriched <- teM_continuous |>
  filter(enrichment_status == "Enriched") |>
  mutate(analysis_group = "Domestication-enriched")

enriched_downregulated <- teM_continuous |>
  filter(enrichment_status == "Enriched", expression_group == "Downregulated") |>
  mutate(analysis_group = "Enriched downregulated")

spearman_one <- function(dat, group_name, carriers_only = FALSE) {
  if (carriers_only) dat <- dat |> filter(any_well_called_pLoF)
  dat <- dat |> filter(is.finite(teM_frequency), is.finite(gene_pLoF_AF))
  if (nrow(dat) < 3) {
    return(tibble(
      analysis_group = group_name,
      subset = if_else(carriers_only, "pLoF carriers", "all genes"),
      n_genes = nrow(dat),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }
  ct <- suppressWarnings(
    cor.test(dat$teM_frequency, dat$gene_pLoF_AF, method = "spearman", exact = FALSE)
  )
  tibble(
    analysis_group = group_name,
    subset = if_else(carriers_only, "pLoF carriers", "all genes"),
    n_genes = nrow(dat),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

correlation_stats <- bind_rows(
  spearman_one(all_teM_positive, "All genes with teM >0"),
  spearman_one(domestication_enriched, "Domestication-enriched"),
  spearman_one(enriched_downregulated, "Enriched downregulated"),
  spearman_one(all_teM_positive, "All genes with teM >0", TRUE),
  spearman_one(domestication_enriched, "Domestication-enriched", TRUE),
  spearman_one(enriched_downregulated, "Enriched downregulated", TRUE)
)

model_presence <- glm(
  any_well_called_pLoF ~ teM_frequency_10pct + log_total_variants + reference_model_warning,
  data = teM_continuous,
  family = binomial()
)

model_af_carriers <- glm(
  gene_pLoF_AF ~ teM_frequency_10pct + log_total_variants + reference_model_warning,
  data = teM_continuous |> filter(any_well_called_pLoF),
  family = quasibinomial(link = "logit")
)

model_results <- bind_rows(
  model_terms(model_presence, "Probability of carrying any well-called pLoF"),
  model_terms(model_af_carriers, "Maximum pLoF AF among pLoF carriers")
)

write_tsv(correlation_stats, file.path(outDir, "suppFig-15b_teM-frequency_pLoF-correlations.tsv"))
write_tsv(model_results, file.path(outDir, "suppFig-15b_teM-frequency_pLoF-models.tsv"))

set.seed(random_seed)
reference_ci <- bootstrap_mean_ci(reference_sample$gene_pLoF_AF)
reference_summary <- tibble(
  n_genes = n_distinct(reference_sample$gene_id),
  n_with_pLoF = sum(reference_sample$any_well_called_pLoF),
  mean_gene_pLoF_AF = reference_ci$mean,
  CI95_low = reference_ci$CI95_low,
  CI95_high = reference_ci$CI95_high
)
write_tsv(reference_summary, file.path(outDir, "suppFig-15b_random-reference_pLoF-AF.tsv"))

presence_term <- model_results |>
  filter(
    model == "Probability of carrying any well-called pLoF",
    term == "teM_frequency_10pct"
  )

carrier_rho <- correlation_stats |>
  filter(
    analysis_group == "All genes with teM >0",
    subset == "pLoF carriers"
  )

subtitle8 <- if (nrow(presence_term) == 1 && nrow(carrier_rho) == 1) {
  paste0(
    "N=", n_distinct(teM_continuous$gene_id),
    ". Per 10% teM-frequency increase, pLoF-presence OR=",
    number(presence_term$exponentiated_estimate, accuracy = 0.01),
    " [", number(presence_term$CI95_low, accuracy = 0.01),
    ", ", number(presence_term$CI95_high, accuracy = 0.01),
    "], P=", format_p(presence_term$p_value),
    "; among pLoF carriers, Spearman ρ=",
    number(carrier_rho$rho, accuracy = 0.01),
    ", P=", format_p(carrier_rho$p_value), "."
  )
} else {
  paste0("Genes observed as teM in at least one accession: N=", n_distinct(teM_continuous$gene_id), ".")
}

p8 <- ggplot() +
  geom_rect(
    data = reference_summary,
    aes(xmin = 0, xmax = 1, ymin = CI95_low, ymax = CI95_high),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.20
  ) +
  geom_hline(
    data = reference_summary,
    aes(yintercept = mean_gene_pLoF_AF),
    inherit.aes = FALSE,
    colour = "grey35",
    linetype = "longdash",
    linewidth = 0.65
  ) +
  geom_point(
    data = all_teM_positive,
    aes(x = teM_frequency, y = gene_pLoF_AF),
    colour = "grey35",
    alpha = 0.075,
    size = 0.85
  ) +
  geom_point(
    data = domestication_enriched,
    aes(x = teM_frequency, y = gene_pLoF_AF),
    colour = "#820a86",
    alpha = 0.23,
    size = 1.15
  ) +
  geom_point(
    data = enriched_downregulated,
    aes(x = teM_frequency, y = gene_pLoF_AF),
    colour = "#d95f5f",
    alpha = 0.72,
    size = 1.55
  ) +
  geom_smooth(
    data = all_teM_positive,
    aes(x = teM_frequency, y = gene_pLoF_AF, colour = "All genes with teM >0"),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 5),
    method.args = list(family = quasibinomial(link = "logit")),
    se = TRUE,
    linewidth = 1.0
  ) +
  geom_smooth(
    data = domestication_enriched,
    aes(x = teM_frequency, y = gene_pLoF_AF, colour = "Domestication-enriched"),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 4),
    method.args = list(family = quasibinomial(link = "logit")),
    se = TRUE,
    linewidth = 1.0
  ) +
  geom_smooth(
    data = enriched_downregulated,
    aes(x = teM_frequency, y = gene_pLoF_AF, colour = "Enriched downregulated"),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 3),
    method.args = list(family = quasibinomial(link = "logit")),
    se = TRUE,
    linewidth = 1.0
  ) +
  scale_colour_manual(values = frequency_colours, drop = FALSE) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = label_percent(accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = label_percent(accuracy = 1)
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Population teM frequency",
    y = "Maximum pLoF allele frequency per gene",
    colour = NULL,
    title = "pLoF frequency versus teM frequency",
    subtitle = subtitle8,
    caption = paste0(
      "The dashed line and grey band show the ",
      "random-reference mean and bootstrap 95% CI."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(1.5, "lines"),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 8, hjust = 0)
  )

save_plot(
  p8,
  "suppFig-15b_pLoF-frequency_vs_teM-frequency.pdf",
  10.5,
  6.8
)
