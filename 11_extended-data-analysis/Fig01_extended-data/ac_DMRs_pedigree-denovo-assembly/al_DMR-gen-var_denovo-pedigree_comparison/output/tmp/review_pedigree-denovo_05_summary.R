library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
root <- "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/ag_ONT-DMR-genetic-vars-comparison"

sv_tabdir    <- file.path(root, "07_sv_state_based", "01_tables")
small_tabdir <- file.path(root, "08_smallvar_explore", "01_tables")

# full DMR universe
dmr_all <- fread(file.path(sv_tabdir, "01_dmr_states.tsv")) %>%
  as_tibble() %>%
  distinct(dmr_id, DMR, change_type, dir_type)

# SV per-DMR summary
sv_dmr <- fread(file.path(sv_tabdir, "06_dmr_sv_dmr_level_summary.tsv")) %>%
  as_tibble() %>%
  select(dmr_id, any_explained_sv = any_explained)

# SNP per-DMR summary
small_dmr <- fread(file.path(small_tabdir, "06_dmr_small_dmr_level_summary.tsv")) %>%
  as_tibble() %>%
  select(dmr_id, any_explained_snp = any_explained)

# merge back to all DMRs
dmr_final <- dmr_all %>%
  left_join(sv_dmr, by = "dmr_id") %>%
  left_join(small_dmr, by = "dmr_id") %>%
  mutate(
    any_explained_sv  = ifelse(is.na(any_explained_sv),  FALSE, any_explained_sv),
    any_explained_snp = ifelse(is.na(any_explained_snp), FALSE, any_explained_snp),
    explain_source = case_when(
      any_explained_sv & any_explained_snp ~ "Both",
      any_explained_sv & !any_explained_snp ~ "SV_only",
      !any_explained_sv & any_explained_snp ~ "SNP_only",
      TRUE ~ "None"
    )
  )

# save full merged table
fwrite(as.data.table(dmr_final),
       file.path(root, "09_dmr_explained_by_sv_snp_or_none.tsv"),
       sep = "\t")

final_counts_full <- dmr_final %>%
  count(DMR, change_type, dir_type, explain_source, name = "n_dmr") %>%
  group_by(DMR, change_type, dir_type) %>%
  mutate(frac = n_dmr / sum(n_dmr)) %>%
  ungroup()

#
plot_df <- final_counts_full %>%
  mutate(
    explain_source = factor(explain_source, levels = c("None", "SV_only", "SNP_only", "Both")),
    dir_type = factor(dir_type, levels = c("loss", "gain", "reversion")),
    change_type = factor(change_type, levels = c("unique", "recurrent"))
  ) %>%
  complete(
    DMR, change_type, dir_type, explain_source,
    fill = list(n_dmr = 0, frac = 0)
  )

p <- ggplot(plot_df, aes(x = dir_type, y = frac, fill = explain_source)) +
  geom_col(width = 0.75) +
  facet_grid(DMR ~ change_type) +
  labs(x = NULL, y = "Fraction of DMRs", fill = NULL) +
  theme_bw(base_size = 12)

ggsave(file.path(root, "00.1_dmr_explained_by_sv_snp_none_by_change_dir.pdf"),
       p, width = 8, height = 5)
#
p_counts <- ggplot(plot_df, aes(x = dir_type, y = n_dmr, fill = explain_source)) +
  geom_col(width = 0.75) +
  facet_grid(DMR ~ change_type, scales = "free_y") +
  labs(x = NULL, y = "Number of DMRs", fill = NULL) +
  theme_bw(base_size = 12)

ggsave(file.path(root, "00.2_dmr_explained_counts_by_sv_snp_none_by_change_dir.pdf"),
       p_counts, width = 8, height = 5)
