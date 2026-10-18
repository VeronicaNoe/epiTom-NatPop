suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# 00) Paths
# =========================================================
root <- "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis"
dmr_file <- "/mnt/disk2/vibanez/10_data-analysis/Fig2/results/02.1_gain-loss_rate_density_per_row.tsv"
small_vcf <- file.path(root, "aa_results", "all-samples.deepvariant.merged.vcf.gz")

outroot <- file.path(root, "ab_DMR-genetic-vars-comparison", "08_smallvar_explore")
tabdir  <- file.path(outroot, "01_tables")
plotdir <- file.path(outroot, "02_plots")

dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

samples <- c("G0-C3", "G0-C6", "G0-C9", "G5-C1", "G5-C2", "G5-C3")

# simple yes/no thresholds
min_support_yes <- 4L
min_support_frac_yes <- 0.75

# =========================================================
# 01) Helpers
# =========================================================
gt_to_bin <- function(x) {
  gt <- sub(":.*$", "", as.character(x))
  gt[gt %in% c(".", "./.", ".|.", "")] <- NA_character_
  out <- rep(NA_integer_, length(gt))
  out[gt %in% c("0/0", "0|0")] <- 0L
  out[!is.na(gt) & !(gt %in% c("0/0", "0|0"))] <- 1L
  out
}

alt_max_nchar <- function(alt) {
  vapply(
    strsplit(as.character(alt), ",", fixed = TRUE),
    function(xx) max(nchar(xx)),
    integer(1)
  )
}

variant_class_from_state <- function(v) {
  n_g0 <- sum(v[1:3] == 1L, na.rm = TRUE)
  n_g5 <- sum(v[4:6] == 1L, na.rm = TRUE)
  
  if (n_g0 > 0 && n_g5 == 0) {
    if (n_g0 == 3) return("fixed_G0")
    return("seg_G0")
  }
  if (n_g5 > 0 && n_g0 == 0) {
    if (n_g5 == 3) return("fixed_G5")
    return("seg_G5")
  }
  if (n_g0 == 3 && n_g5 == 3) return("shared_all")
  if (n_g0 > 0 && n_g5 > 0) return("mixed_both")
  return("none")
}

state_signature <- function(dt, cols) {
  apply(dt[, ..cols], 1, function(z) {
    z <- ifelse(is.na(z), "N", as.character(z))
    paste(z, collapse = "")
  })
}

save_plot <- function(p, prefix, width = 7, height = 5) {
  ggsave(file.path(plotdir, paste0(prefix, ".pdf")), plot = p, width = width, height = height)
  ggsave(file.path(plotdir, paste0(prefix, ".png")), plot = p, width = width, height = height, dpi = 300)
}

# =========================================================
# 02) DMR states
# =========================================================
dmr <- fread(dmr_file, sep = "\t", header = TRUE)

dmr_states <- dmr[, .(
  dmr_id = paste(chr, start, end, DMR, sep = ":"),
  dmr_chrom = chr,
  dmr_start = as.integer(start),
  dmr_end   = as.integer(end),
  dmr_len   = as.integer(end - start + 1L),
  DMR,
  change_type,
  dir_type,
  ancestral_state,
  derived_state,
  `dmr_G0-C3` = as.integer(`bin_C3-G0`),
  `dmr_G0-C6` = as.integer(`bin_C6-G0`),
  `dmr_G0-C9` = as.integer(`bin_C9-G0`),
  `dmr_G5-C1` = as.integer(`bin_C1-G5`),
  `dmr_G5-C2` = as.integer(`bin_C2-G5`),
  `dmr_G5-C3` = as.integer(`bin_C3-G5`)
)]

fwrite(dmr_states, file.path(tabdir, "01_dmr_states.tsv"), sep = "\t")

# =========================================================
# 03) Small-variant states
# =========================================================
small <- fread(
  cmd = sprintf("zcat %s | grep -v '^##'", shQuote(small_vcf)),
  sep = "\t",
  header = TRUE
)

small_states <- small[, .(
  chrom = `#CHROM`,
  pos   = as.integer(POS),
  ref   = REF,
  alt   = ALT
)]

for (s in samples) {
  small_states[, (paste0("gt_", s)) := sub(":.*$", "", small[[s]])]
}
for (s in samples) {
  small_states[, (paste0("var_", s)) := gt_to_bin(get(paste0("gt_", s)))]
}

small_states[, start := pos]
small_states[, end := as.integer(pos + pmax(nchar(ref), alt_max_nchar(alt)) - 1L)]
small_states[, var_len := as.integer(end - start + 1L)]
small_states[, variant_id := paste(chrom, pos, ref, alt, sep = ":")]
small_states[, variant_type := ifelse(
  nchar(ref) == 1L & alt_max_nchar(alt) == 1L & !grepl(",", alt, fixed = TRUE),
  "SNP", "INDEL"
)]

small_states[, n_g0 := rowSums(.SD == 1L, na.rm = TRUE), .SDcols = c("var_G0-C3","var_G0-C6","var_G0-C9")]
small_states[, n_g5 := rowSums(.SD == 1L, na.rm = TRUE), .SDcols = c("var_G5-C1","var_G5-C2","var_G5-C3")]
small_states[, variant_class := apply(.SD, 1, variant_class_from_state),
             .SDcols = c("var_G0-C3","var_G0-C6","var_G0-C9","var_G5-C1","var_G5-C2","var_G5-C3")]

fwrite(small_states, file.path(tabdir, "02_small_states.tsv"), sep = "\t")

# =========================================================
# 04) Overlap DMRs with small variants
# =========================================================
dmr_int <- copy(dmr_states)[, .(
  chrom = dmr_chrom,
  start = dmr_start,
  end   = dmr_end,
  
  dmr_id, dmr_chrom, dmr_start, dmr_end, dmr_len,
  DMR, change_type, dir_type, ancestral_state, derived_state,
  `dmr_G0-C3`, `dmr_G0-C6`, `dmr_G0-C9`,
  `dmr_G5-C1`, `dmr_G5-C2`, `dmr_G5-C3`
)]
setkey(dmr_int, chrom, start, end)

small_int <- copy(small_states)[, .(
  chrom = chrom,
  start = start,
  end   = end,
  
  variant_id,
  variant_chrom = chrom,
  variant_start = start,
  variant_end   = end,
  variant_pos   = pos,
  variant_ref   = ref,
  variant_alt   = alt,
  variant_type,
  var_len,
  n_g0,
  n_g5,
  variant_class,
  
  `gt_G0-C3`, `gt_G0-C6`, `gt_G0-C9`,
  `gt_G5-C1`, `gt_G5-C2`, `gt_G5-C3`,
  
  `var_G0-C3`, `var_G0-C6`, `var_G0-C9`,
  `var_G5-C1`, `var_G5-C2`, `var_G5-C3`
)]
setkey(small_int, chrom, start, end)

ov_small <- foverlaps(
  x = small_int,
  y = dmr_int,
  type = "any",
  nomatch = 0L
)

ov_small <- ov_small[, .(
  dmr_id, dmr_chrom, dmr_start, dmr_end, dmr_len,
  DMR, change_type, dir_type, ancestral_state, derived_state,
  
  `dmr_G0-C3`, `dmr_G0-C6`, `dmr_G0-C9`,
  `dmr_G5-C1`, `dmr_G5-C2`, `dmr_G5-C3`,
  
  variant_id, variant_chrom, variant_pos, variant_start, variant_end,
  variant_ref, variant_alt, variant_type, var_len,
  n_g0, n_g5, variant_class,
  
  `gt_G0-C3`, `gt_G0-C6`, `gt_G0-C9`,
  `gt_G5-C1`, `gt_G5-C2`, `gt_G5-C3`,
  
  `var_G0-C3`, `var_G0-C6`, `var_G0-C9`,
  `var_G5-C1`, `var_G5-C2`, `var_G5-C3`
)]

ov_small[, overlap_bp := pmin(dmr_end, variant_end) - pmax(dmr_start, variant_start) + 1L]

# =========================================================
# 05) Score DMR–smallvar support
#     DMR NA -> 0 ; only called variant samples count
# =========================================================
dmr_cols <- paste0("dmr_", samples)
var_cols <- paste0("var_", samples)

dmr_m <- as.matrix(ov_small[, ..dmr_cols])
var_m <- as.matrix(ov_small[, ..var_cols])

storage.mode(dmr_m) <- "integer"
storage.mode(var_m) <- "integer"

dmr_m0 <- dmr_m
dmr_m0[is.na(dmr_m0)] <- 0L

var_called <- !is.na(var_m)
n_var_called <- rowSums(var_called)

support_direct <- var_called & (dmr_m0 == var_m)
support_inverse <- var_called & (dmr_m0 != var_m)

n_support_direct <- rowSums(support_direct)
n_support_inverse <- rowSums(support_inverse)

support_direction <- ifelse(n_support_direct >= n_support_inverse, "direct", "inverse")
n_support <- pmax(n_support_direct, n_support_inverse)
n_contradict <- n_var_called - n_support
support_frac <- ifelse(n_var_called > 0, n_support / n_var_called, NA_real_)

ov_small[, n_dmr_na := rowSums(is.na(dmr_m))]
ov_small[, n_var_called := n_var_called]
ov_small[, n_support_direct := n_support_direct]
ov_small[, n_support_inverse := n_support_inverse]
ov_small[, n_support := n_support]
ov_small[, n_contradict := n_contradict]
ov_small[, support_direction := support_direction]
ov_small[, support_frac := support_frac]

# =========================================================
# 06) Simple explanatory yes/no
# =========================================================
ov_small[, explanatory_yesno := ifelse(
  !is.na(support_frac) &
    n_support >= min_support_yes &
    support_frac >= min_support_frac_yes,
  "yes", "no"
)]

ov_small[, variant_signature := state_signature(.SD, var_cols)]

# =========================================================
# 07) Export pair table
# =========================================================
small_pairs <- ov_small[, .(
  dmr_id, dmr_chrom, dmr_start, dmr_end, dmr_len,
  DMR, change_type, dir_type, ancestral_state, derived_state,
  
  variant_id, variant_chrom, variant_pos, variant_start, variant_end,
  variant_ref, variant_alt, variant_type, var_len,
  n_g0, n_g5, variant_class, variant_signature,
  
  overlap_bp,
  
  `dmr_G0-C3`, `dmr_G0-C6`, `dmr_G0-C9`,
  `dmr_G5-C1`, `dmr_G5-C2`, `dmr_G5-C3`,
  
  `gt_G0-C3`, `gt_G0-C6`, `gt_G0-C9`,
  `gt_G5-C1`, `gt_G5-C2`, `gt_G5-C3`,
  
  `var_G0-C3`, `var_G0-C6`, `var_G0-C9`,
  `var_G5-C1`, `var_G5-C2`, `var_G5-C3`,
  
  n_dmr_na,
  n_var_called,
  n_support_direct,
  n_support_inverse,
  n_support,
  n_contradict,
  support_direction,
  support_frac,
  explanatory_yesno
)]

fwrite(small_pairs, file.path(tabdir, "03_dmr_small_overlap_pairs.tsv"), sep = "\t")

# =========================================================
# 08) Condense by DMR + signature
# =========================================================
small_sig <- small_pairs %>%
  as_tibble() %>%
  group_by(
    dmr_id, DMR, change_type, dir_type,
    variant_signature, variant_class, variant_type,
    support_direction, explanatory_yesno
  ) %>%
  summarise(
    n_var = n(),
    max_support = max(n_support, na.rm = TRUE),
    max_support_frac = max(support_frac, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(as.data.table(small_sig), file.path(tabdir, "04_dmr_small_signature_summary.tsv"), sep = "\t")

# =========================================================
# 09) Best small variant per DMR
# =========================================================
small_best <- copy(small_pairs)

small_best[, explanatory_rank := ifelse(explanatory_yesno == "yes", 1L, 0L)]

setorder(
  small_best,
  dmr_id,
  -explanatory_rank,
  -n_support,
  -support_frac,
  n_contradict,
  variant_id
)

small_best <- small_best[, .SD[1], by = dmr_id]
small_best[, explanatory_rank := NULL]

fwrite(small_best, file.path(tabdir, "05_dmr_small_best_per_dmr.tsv"), sep = "\t")

# =========================================================
# 10) DMR-level summary
# =========================================================
small_dmr_summary <- small_pairs %>%
  as_tibble() %>%
  group_by(dmr_id, DMR, change_type, dir_type) %>%
  summarise(
    n_overlap_small = n(),
    n_explanatory_small = sum(explanatory_yesno == "yes", na.rm = TRUE),
    any_explained = n_explanatory_small > 0,
    n_explanatory_signatures = n_distinct(variant_signature[explanatory_yesno == "yes"]),
    .groups = "drop"
  ) %>%
  left_join(
    small_best %>%
      as_tibble() %>%
      select(
        dmr_id,
        best_variant_id = variant_id,
        best_variant_type = variant_type,
        best_variant_class = variant_class,
        best_n_support = n_support,
        best_support_frac = support_frac,
        best_support_direction = support_direction,
        best_explanatory = explanatory_yesno
      ),
    by = "dmr_id"
  )

fwrite(as.data.table(small_dmr_summary), file.path(tabdir, "06_dmr_small_dmr_level_summary.tsv"), sep = "\t")

# =========================================================
# 11) Summary counts
# =========================================================
small_counts <- small_dmr_summary %>%
  count(DMR, change_type, dir_type, any_explained, name = "n_dmr") %>%
  group_by(DMR, change_type, dir_type) %>%
  mutate(frac = n_dmr / sum(n_dmr)) %>%
  ungroup()

small_best_class_counts <- small_dmr_summary %>%
  filter(best_explanatory == "yes") %>%
  count(DMR, change_type, dir_type, best_variant_class, name = "n_dmr") %>%
  arrange(DMR, change_type, dir_type, best_variant_class)

fwrite(as.data.table(small_counts), file.path(tabdir, "07_small_explained_counts_binary.tsv"), sep = "\t")
fwrite(as.data.table(small_best_class_counts), file.path(tabdir, "08_small_best_explanatory_class_counts.tsv"), sep = "\t")

# =========================================================
# 12) Plots
# =========================================================
plot_df <- small_counts %>%
  mutate(any_explained = ifelse(any_explained, "Explained", "Not explained"))

p1 <- ggplot(plot_df, aes(x = dir_type, y = frac, fill = any_explained)) +
  geom_col(width = 0.75) +
  facet_grid(DMR ~ change_type) +
  labs(x = NULL, y = "Fraction of DMRs", fill = NULL) +
  theme_bw(base_size = 12)

save_plot(p1, "01_small_explained_fraction_binary", width = 8, height = 5)

p2 <- ggplot(small_best_class_counts, aes(x = best_variant_class, y = n_dmr)) +
  geom_col() +
  facet_grid(DMR ~ change_type + dir_type, scales = "free_y") +
  labs(x = "Best explanatory small-variant class", y = "Number of DMRs") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p2, "02_small_best_explanatory_class_counts", width = 12, height = 6)

p3 <- ggplot(small_dmr_summary, aes(x = n_overlap_small)) +
  geom_histogram(binwidth = 1) +
  facet_grid(DMR ~ change_type) +
  labs(x = "Number of overlapping small variants per DMR", y = "Number of DMRs") +
  theme_bw(base_size = 11)

save_plot(p3, "03_small_n_overlap_histogram", width = 8, height = 5)

cat("Done.\n")
cat("Tables saved in: ", tabdir, "\n")
cat("Plots saved in:  ", plotdir, "\n")