suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# 00) Paths and settings
# =========================================================
root <- "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis"
dmr_file <- "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/af_DMR-epi-rate/03.2_gain-loss_rate_density_per_row.tsv"
sv_vcf   <- file.path(root, "aa_results", "cervil-SV-ONT.PASS.SUP5.tags.vcf.gz")
small_vcf <- file.path(root, "aa_results", "cervil-SNPs-ONT.GQ20_DP5.tags.vcf.")
outroot  <- file.path(root, "ag_ONT-DMR-genetic-vars-comparison", "07_sv_state_based")
tabdir   <- file.path(outroot, "01_tables")
plotdir  <- file.path(outroot, "02_plots")

dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

samples <- c("G0-C3", "G0-C6", "G0-C9", "G5-C1", "G5-C2", "G5-C3")

# overlap thresholds
min_overlap_frac_strong   <- 0.80
min_overlap_frac_possible <- 0.50

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

extract_info_value <- function(info, key) {
  sapply(strsplit(info, ";", fixed = TRUE), function(xx) {
    hit <- xx[startsWith(xx, paste0(key, "="))]
    if (length(hit) == 0) return(NA_character_)
    sub(paste0("^", key, "="), "", hit[1])
  }, USE.NAMES = FALSE)
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

save_plot <- function(p, prefix, width = 7, height = 5) {
  ggsave(file.path(plotdir, paste0(prefix, ".pdf")), plot = p, width = width, height = height)
  ggsave(file.path(plotdir, paste0(prefix, ".png")), plot = p, width = width, height = height, dpi = 300)
}

alt_max_nchar <- function(alt) {
  vapply(
    strsplit(as.character(alt), ",", fixed = TRUE),
    function(xx) max(nchar(xx)),
    integer(1)
  )
}


# =========================================================
# 02) Read DMR table and build DMR states
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
# 03) Read SV VCF and build SV states
# =========================================================
sv <- fread(cmd = sprintf("zcat %s | grep -v '^##'", shQuote(sv_vcf)),sep = "\t",header = TRUE)

sv_states <- sv[, .(chrom = `#CHROM`,pos   = as.integer(POS),ref   = REF, alt   = ALT,
  info  = INFO)]

sv_states[, variant_end  := as.integer(extract_info_value(info, "END"))]
sv_states[, variant_type := extract_info_value(info, "SVTYPE")]
sv_states[is.na(variant_end), variant_end := pos]

# keep raw GTs too, for inspection
for (s in samples) {
  sv_states[, (paste0("gt_", s)) := sub(":.*$", "", sv[[s]])]
}
for (s in samples) {
  sv_states[, (paste0("var_", s)) := gt_to_bin(get(paste0("gt_", s)))]
}

sv_states[, start := pos]
sv_states[, end   := variant_end]
sv_states[, sv_len := as.integer(end - start + 1L)]
sv_states[, variant_id := paste(chrom, pos, variant_type, variant_end, sep = ":")]

sv_states[, n_g0 := rowSums(.SD == 1L, na.rm = TRUE), .SDcols = c("var_G0-C3","var_G0-C6","var_G0-C9")]
sv_states[, n_g5 := rowSums(.SD == 1L, na.rm = TRUE), .SDcols = c("var_G5-C1","var_G5-C2","var_G5-C3")]
sv_states[, variant_class := apply(.SD, 1, variant_class_from_state),
          .SDcols = c("var_G0-C3","var_G0-C6","var_G0-C9","var_G5-C1","var_G5-C2","var_G5-C3")]

fwrite(sv_states, file.path(tabdir, "02_sv_states.tsv"), sep = "\t")

# =========================================================
# 04) Overlap DMRs with SVs
# =========================================================
dmr_int <- copy(dmr_states)[, .(chrom = dmr_chrom,start = dmr_start,end   = dmr_end,
                                dmr_id, dmr_chrom, dmr_start, dmr_end, dmr_len,
                                DMR, change_type, dir_type, ancestral_state, derived_state,
                                `dmr_G0-C3`, `dmr_G0-C6`, `dmr_G0-C9`,
                                `dmr_G5-C1`, `dmr_G5-C2`, `dmr_G5-C3`)]
setkey(dmr_int, chrom, start, end)

sv_int <- copy(sv_states)[, .(chrom = chrom,start = start,end   = end,
                              variant_id, variant_chrom = chrom,
                              variant_start = start, variant_end   = end,
                              variant_pos   = pos,  variant_ref   = ref,
                              variant_alt   = alt,  variant_type,  sv_len,
                              n_g0,  n_g5,  variant_class,  `gt_G0-C3`, 
                              `gt_G0-C6`, `gt_G0-C9`,  `gt_G5-C1`, `gt_G5-C2`, 
                              `gt_G5-C3`,  `var_G0-C3`, `var_G0-C6`, `var_G0-C9`,
                              `var_G5-C1`, `var_G5-C2`, `var_G5-C3`)]
setkey(sv_int, chrom, start, end)

ov_sv <- foverlaps( x = sv_int, y = dmr_int, type = "any", nomatch = 0L)

ov_sv <- ov_sv[, .(dmr_id, dmr_chrom, dmr_start, dmr_end, dmr_len,
                   DMR, change_type, dir_type, ancestral_state, derived_state,
                   `dmr_G0-C3`, `dmr_G0-C6`, `dmr_G0-C9`,
                   `dmr_G5-C1`, `dmr_G5-C2`, `dmr_G5-C3`,
                   variant_id, variant_chrom, variant_pos, variant_start, variant_end,
                   variant_ref, variant_alt, variant_type, sv_len,
                   n_g0, n_g5, variant_class,  `gt_G0-C3`, `gt_G0-C6`, `gt_G0-C9`,
                   `gt_G5-C1`, `gt_G5-C2`, `gt_G5-C3`,  `var_G0-C3`, `var_G0-C6`, `var_G0-C9`,
                   `var_G5-C1`, `var_G5-C2`, `var_G5-C3`)]

ov_sv[, overlap_bp   := pmin(dmr_end, variant_end) - pmax(dmr_start, variant_start) + 1L]
ov_sv[, overlap_frac := overlap_bp / dmr_len]

# =========================================================
# 05) Score DMR–SV explanatory support
# =========================================================
# =========================================================
# 05) Score DMR–SV support
#     DMR NA is treated as 0
#     only samples with called SV genotype are counted
# =========================================================
dmr_cols <- paste0("dmr_", samples)
var_cols <- paste0("var_", samples)

dmr_m <- as.matrix(ov_sv[, ..dmr_cols])
var_m <- as.matrix(ov_sv[, ..var_cols])

storage.mode(dmr_m) <- "integer"
storage.mode(var_m) <- "integer"

# treat DMR NA as 0
dmr_m0 <- dmr_m
dmr_m0[is.na(dmr_m0)] <- 0L

# only SV-called samples contribute
sv_called <- !is.na(var_m)
n_sv_called <- rowSums(sv_called)

# support under direct and inverse models
support_direct <- sv_called & (dmr_m0 == var_m)
support_inverse <- sv_called & (dmr_m0 != var_m)

n_support_direct <- rowSums(support_direct)
n_support_inverse <- rowSums(support_inverse)

support_direction <- ifelse(n_support_direct >= n_support_inverse, "direct", "inverse")
n_support <- pmax(n_support_direct, n_support_inverse)
n_contradict <- n_sv_called - n_support
support_frac <- ifelse(n_sv_called > 0, n_support / n_sv_called, NA_real_)

# annotation only
n_dmr_na <- rowSums(is.na(dmr_m))

ov_sv[, n_sv_called := n_sv_called]
ov_sv[, n_support_direct := n_support_direct]
ov_sv[, n_support_inverse := n_support_inverse]
ov_sv[, n_support := n_support]
ov_sv[, n_contradict := n_contradict]
ov_sv[, support_direction := support_direction]
ov_sv[, support_frac := support_frac]
ov_sv[, n_dmr_na := n_dmr_na]

# =========================================================
# 06) Simple explanatory yes/no
#     support-based, class kept as annotation
# =========================================================
min_support_yes <- 4L
min_support_frac_yes <- 0.75

ov_sv[, explanatory_yesno := ifelse(
  !is.na(support_frac) &
    n_support >= min_support_yes &
    support_frac >= min_support_frac_yes,
  "yes", "no"
)]

# =========================================================
# 07) Export pair table
# =========================================================
sv_pairs <- ov_sv[, .(
  dmr_id, dmr_chrom, dmr_start, dmr_end, dmr_len,
  DMR, change_type, dir_type, ancestral_state, derived_state,
  
  variant_id, variant_chrom, variant_pos, variant_start, variant_end,
  variant_ref, variant_alt, variant_type, sv_len,
  n_g0, n_g5, variant_class,
  
  overlap_bp, overlap_frac,
  
  `dmr_G0-C3`, `dmr_G0-C6`, `dmr_G0-C9`,
  `dmr_G5-C1`, `dmr_G5-C2`, `dmr_G5-C3`,
  
  `gt_G0-C3`, `gt_G0-C6`, `gt_G0-C9`,
  `gt_G5-C1`, `gt_G5-C2`, `gt_G5-C3`,
  
  `var_G0-C3`, `var_G0-C6`, `var_G0-C9`,
  `var_G5-C1`, `var_G5-C2`, `var_G5-C3`,
  
  n_dmr_na,
  n_sv_called,
  n_support_direct,
  n_support_inverse,
  n_support,
  n_contradict,
  support_direction,
  support_frac,
  explanatory_yesno
)]

fwrite(sv_pairs, file.path(tabdir, "03_dmr_sv_overlap_pairs.tsv"), sep = "\t")

# one row per DMR–SV–sample
sv_long <- rbindlist(lapply(samples, function(s) {
  sv_pairs[, .(
    dmr_id, dmr_chrom, dmr_start, dmr_end, DMR, change_type, dir_type,
    variant_id, variant_type, variant_class,
    sample = s,
    dmr_state = get(paste0("dmr_", s)),
    gt_state  = get(paste0("gt_", s)),
    sv_state  = get(paste0("var_", s)),
    support_direction, n_support, support_frac, explanatory_yesno
  )]
}), use.names = TRUE)

fwrite(sv_long, file.path(tabdir, "04_dmr_sv_overlap_long.tsv"), sep = "\t")

# =========================================================
# 08) Best SV per DMR
# =========================================================
sv_best <- copy(sv_pairs)

sv_best[, explanatory_rank := ifelse(explanatory_yesno == "yes", 1L, 0L)]

setorder(
  sv_best,
  dmr_id,
  -explanatory_rank,
  -n_support,
  -support_frac,
  n_contradict,
  -overlap_frac,
  sv_len,
  variant_id
)

sv_best <- sv_best[, .SD[1], by = dmr_id]
sv_best[, explanatory_rank := NULL]

fwrite(sv_best, file.path(tabdir, "05_dmr_sv_best_per_dmr.tsv"), sep = "\t")

# =========================================================
# 09) DMR-level summary
# =========================================================
sv_dmr_summary <- sv_pairs %>%
  as_tibble() %>%
  group_by(dmr_id, DMR, change_type, dir_type) %>%
  summarise(
    n_overlap_sv = n(),
    n_explanatory_sv = sum(explanatory_yesno == "yes", na.rm = TRUE),
    any_explained = n_explanatory_sv > 0,
    .groups = "drop"
  ) %>%
  left_join(
    sv_best %>%
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

fwrite(as.data.table(sv_dmr_summary), file.path(tabdir, "06_dmr_sv_dmr_level_summary.tsv"), sep = "\t")

# =========================================================
# 10) Summary counts
# =========================================================
sv_counts <- sv_dmr_summary %>%
  count(DMR, change_type, dir_type, any_explained, name = "n_dmr") %>%
  group_by(DMR, change_type, dir_type) %>%
  mutate(frac = n_dmr / sum(n_dmr)) %>%
  ungroup()

sv_best_class_counts <- sv_dmr_summary %>%
  filter(best_explanatory == "yes") %>%
  count(DMR, change_type, dir_type, best_variant_class, name = "n_dmr") %>%
  arrange(DMR, change_type, dir_type, best_variant_class)

fwrite(as.data.table(sv_counts), file.path(tabdir, "07_sv_explained_counts_binary.tsv"), sep = "\t")
fwrite(as.data.table(sv_best_class_counts), file.path(tabdir, "08_sv_best_explanatory_class_counts.tsv"), sep = "\t")

# =========================================================
# 11) Plots
# =========================================================
plot_df <- sv_counts %>%
  mutate(any_explained = ifelse(any_explained, "Explained", "Not explained"))

p1 <- ggplot(plot_df, aes(x = dir_type, y = frac, fill = any_explained)) +
  geom_col(width = 0.75) +
  facet_grid(DMR ~ change_type) +
  labs(x = NULL, y = "Fraction of DMRs", fill = NULL) +
  theme_bw(base_size = 12)

save_plot(p1, "01_sv_explained_fraction_binary", width = 8, height = 5)

p2 <- ggplot(sv_best_class_counts, aes(x = best_variant_class, y = n_dmr)) +
  geom_col() +
  facet_grid(DMR ~ change_type + dir_type, scales = "free_y") +
  labs(x = "Best explanatory SV class", y = "Number of DMRs") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p2, "02_sv_best_explanatory_class_counts", width = 12, height = 6)

cat("Done.\n")
cat("Tables saved in: ", tabdir, "\n")
cat("Plots saved in:  ", plotdir, "\n")
