#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

# =========================================================
# 00) Paths
# =========================================================
ONT_data_dir <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ai_DMR-classification/bd_merged-methylation/"
ONT_data_loci <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ai_DMR-classification/bb_get-collapsed-loci/"

Biseq_data_dir <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation/cervil-denovo-3x/"
Biseq_data_loci <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree/cervil-denovo-3x/"

Biseq_sl25_denovo_dir <- "/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation/cervil-sl25/"

outDir <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/al_DMR-gen-var_denovo-pedigree_comparison/"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

samples <- c("G0-C3", "G0-C6", "G0-C9", "G5-C1", "G5-C2", "G5-C3")
dmr_types <- c("C-DMR", "CG-DMR")
HEATMAP_N <- 5000

dmr_colors <- c("C-DMR" = "#820a86",  "CG-DMR" = "#ffca7b")

dataset_colors <- c(  "ONT" = "#377eb8",  "BIseq" = "#4daf4a")

# =========================================================
# 01) Read lifted SL2.5 -> polished BIseq DMRs
#      Used only to filter ONT
# =========================================================
read_lifted_dmr <- function(dmr_type, lifted_dir) {
  f <- file.path(lifted_dir, paste0(dmr_type, ".SL25_to_polished.bed"))
  
  if (!file.exists(f)) stop("Missing lifted DMR file: ", f)
  
  dt <- fread(f, header = FALSE, sep = "\t", fill = TRUE, data.table = TRUE)
  
  dt <- dt[, .(
    DMR = dmr_type,
    chr = as.character(V1),
    start = as.integer(V2),
    end = as.integer(V3),
    dmr_id = as.character(V4)
  )]
  
  unique(dt)
}

lifted_dt <- rbindlist(lapply(dmr_types, read_lifted_dmr, lifted_dir = Biseq_sl25_denovo_dir),fill = TRUE)
#fwrite(lifted_dt,  file = file.path(outDir, "02.0_lifted_SL25_to_polished_DMRs_used_to_filter_ONT.tsv"), sep = "\t",  quote = FALSE)

# =========================================================
# 02) Helper functions
# =========================================================
filter_to_lifted <- function(dt, lifted_dt) {
  x <- copy(dt)
  y <- copy(lifted_dt)
  if (nrow(x) == 0) return(x)
  x[, row_id := .I]
  setkey(y, DMR, chr, start, end)
  setkey(x, DMR, chr, start, end)
  hits <- foverlaps(
    x,
    y,
    by.x = c("DMR", "chr", "start", "end"),
    by.y = c("DMR", "chr", "start", "end"),
    type = "any",
    nomatch = 0L
  )
  
  if (nrow(hits) == 0) {
    x[, row_id := NULL]
    x[, dmr_id := NA_character_]
    return(x[0])
  }
  
  hit_ids <- hits[, .(
    dmr_id = paste(unique(dmr_id), collapse = ";")
  ), by = row_id]
  
  out <- merge(x, hit_ids, by = "row_id", all = FALSE)
  out[, row_id := NULL]
  
  out
}

calc_density_threshold <- function(row_values) {
  row_values <- row_values[!is.na(row_values)]
  
  if (length(row_values) < 2) return(NA_real_)
  if (length(unique(row_values)) == 1) return(unique(row_values))
  if (length(row_values) == 2) return(mean(row_values))
  
  density_data <- tryCatch(density(row_values), error = function(e) NULL)
  if (is.null(density_data)) return(median(row_values))
  
  dens_diff <- diff(density_data$y)
  valley_points <- which(diff(sign(dens_diff)) == 2)
  
  if (length(valley_points) == 0) return(median(row_values))
  
  density_data$x[valley_points][which.min(density_data$y[valley_points])]
}

safe_range <- function(z) {
  if (all(is.na(z))) return(NA_real_)
  max(z, na.rm = TRUE) - min(z, na.rm = TRUE)
}

# =========================================================
# 03) Read one dataset
#      ONT is filtered to lifted SL2.5->polished DMRs
# =========================================================
read_dmr_dataset <- function(data_dir, loci_dir, dataset_name, filter_with_lifted = FALSE) {
  input <- list.files(
    path = data_dir,
    pattern = "\\.bed$|\\.bed\\.gz$",
    full.names = FALSE  )
  input <- input[grepl("^(C-DMR|CG-DMR)_all.*\\.bed(\\.gz)?$", input)]
  if (length(input) == 0) {
    stop("No merged DMR BED files found in: ", data_dir)
  }
  out_list <- list()
  for (file_name in input) {
    message("Processing ", dataset_name, ": ", file_name)
    dmr <- strsplit(basename(file_name), "_", fixed = TRUE)[[1]][1]
    if (!dmr %in% dmr_types) next
    header_file <- file.path(data_dir, paste0(dmr, "_all.colNames.tsv"))
    loci_file <- file.path(loci_dir, paste0(dmr, "-loci_collapsed.bed"))
    data_file <- file.path(data_dir, file_name)
    if (!file.exists(header_file)) stop("Missing header file: ", header_file)
    header <- fread(
      header_file,
      sep = "\t",
      header = FALSE,
      data.table = TRUE,
      fill = TRUE,
      na.strings = c(""),
      nThread = 20
    )
    tmp <- fread(
      data_file,
      sep = "\t",
      header = FALSE,
      data.table = TRUE,
      fill = TRUE,
      na.strings = c(""),
      nThread = 20
    )
    new_names <- c("chr", "start", "end", header$V1)
    if (length(new_names) != ncol(tmp)) {
      stop(
        "Column number mismatch for ", data_file,
        "\nExpected: ", length(new_names),
        "\nObserved: ", ncol(tmp)
      )
    }
    setnames(tmp, new_names)
    valid_cols <- intersect(samples, colnames(tmp))
    tmp <- tmp[, c("chr", "start", "end", valid_cols), with = FALSE]
    tmp[, `:=`(
      chr = as.character(chr),
      start = as.integer(start),
      end = as.integer(end),
      DMR = dmr,
      dataset = dataset_name,
      source_file = file_name
    )]
    tmp[, (valid_cols) := lapply(.SD, as.numeric), .SDcols = valid_cols]
    # Keep only collapsed loci if available
    if (file.exists(loci_file)) {
      dmr_loci <- fread(
        loci_file,
        sep = "\t",
        header = FALSE,
        data.table = TRUE,
        fill = TRUE,
        na.strings = c(""),
        nThread = 20
      )
      setnames(dmr_loci, old = c("V1", "V2"), new = c("chr", "start"))
      dmr_loci <- unique(dmr_loci[, .(
        chr = as.character(chr),
        start = as.integer(start)
      )])
      tmp <- merge(tmp, dmr_loci, by = c("chr", "start"), all = FALSE)
    }
    # Filter only ONT with lifted SL2.5 
    if (filter_with_lifted) {
      tmp <- filter_to_lifted(tmp, lifted_dt[DMR == dmr])
    } else {
      tmp[, dmr_id := NA_character_]
    }
    if (nrow(tmp) == 0) {
      message("No rows retained for ", dataset_name, " ", dmr)
      next
    }
    tmp[, locus_id := fifelse(
      !is.na(dmr_id) & dmr_id != "",
      dmr_id,
      paste(chr, DMR, start, sep = "_")
    )]
    meth_mat <- as.matrix(tmp[, ..valid_cols])
    tmp[, density_threshold := apply(meth_mat, 1, calc_density_threshold)]
    for (col in valid_cols) {
      tmp[, paste0("bin_", col) := fifelse(
        is.na(get(col)) | is.na(density_threshold),
        NA_integer_,
        fifelse(get(col) > density_threshold, 1L, 0L)
      )]
    }
    bin_cols <- paste0("bin_", valid_cols)
    tmp[, `:=`(
      count_1 = rowSums(.SD == 1, na.rm = TRUE),
      count_0 = rowSums(.SD == 0, na.rm = TRUE),
      count_NA = rowSums(is.na(.SD))
    ), .SDcols = bin_cols]
    
    tmp[, ancestral_state := fifelse(
      count_1 > count_0,
      1L,
      fifelse(count_0 > count_1, 0L, 2L)
    )]
    
    tmp[, derived_state := fifelse(
      ancestral_state == 1L,
      0L,
      fifelse(ancestral_state == 0L, 1L, 2L)
    )]
    
    tmp[, change_type := fifelse(
      (count_0 == 1L & count_1 > 1L) | (count_0 > 1L & count_1 == 1L),
      "unique",
      "recurrent"
    )]
    
    tmp[, dir_type := fifelse(
      derived_state == 1L,
      "gain",
      fifelse(derived_state == 0L, "loss", "reversion")
    )]
    
    out_list[[length(out_list) + 1]] <- tmp
  }
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

df_ont <- read_dmr_dataset(data_dir = ONT_data_dir, loci_dir = ONT_data_loci, dataset_name = "ONT",  filter_with_lifted = TRUE)

df_biseq <- read_dmr_dataset(data_dir = Biseq_data_dir, loci_dir = Biseq_data_loci, dataset_name = "BIseq",filter_with_lifted = FALSE)

df <- rbindlist(list(df_ont, df_biseq), use.names = TRUE, fill = TRUE)
setorder(df, dataset, DMR, chr, start, end)

fwrite(df, file = file.path(outDir, "02.1_ONT_BIseq_DMR_methylation_binary_states.tsv"),sep = "\t",  quote = FALSE)

# =========================================================
# 04) Gain/loss/reversion summary, both datasets
# =========================================================
filtered_data <- df[count_NA == 0]

dmr_summary <- filtered_data[, .N, by = .(dataset, DMR, change_type, dir_type)]
dmr_summary[, proportion := N / sum(N), by = .(dataset, DMR)]
dmr_summary[, change_dir := paste(change_type, dir_type, sep = ":")]

desired_order <- c(
  "unique:loss",
  "unique:gain",
  "recurrent:loss",
  "recurrent:gain",
  "recurrent:reversion"
)

dmr_summary[, change_dir := factor(change_dir, levels = desired_order)]

fwrite(dmr_summary,file = file.path(outDir, "02.3_gain_loss_reversion_summary_ONT_filtered_BIseq_unfiltered.tsv"),sep = "\t",  quote = FALSE)

p_summary <- ggplot(dmr_summary, aes(x = change_dir, y = proportion, fill = DMR, color = DMR)) +
  geom_col(position = "dodge") +
  facet_grid(DMR ~ dataset, scales = "free_y") +
  scale_color_manual(values = dmr_colors) +
  scale_fill_manual(values = dmr_colors) +
  labs(title = "DMR state changes",
    x = "Change type:direction",
    y = "Proportion"  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none"  )

ggsave(file.path(outDir, "02.3_gain_loss_reversion_summary_ONT_BIseq.pdf"),p_summary,  width = 10,  height = 6)

# =========================================================
# 06) Synthetic methylation-state VCFs, one per dataset and DMR
# =========================================================
bin_to_gt <- function(x) {
  out <- rep("./.", length(x))
  out[!is.na(x) & x == 0] <- "0/0"
  out[!is.na(x) & x == 1] <- "1/1"
  out
}

write_mvcf <- function(dt, dataset_name, dmr_type, out_prefix, bgzip = TRUE) {
  x <- copy(dt[dataset == dataset_name & DMR == dmr_type])
  
  if (nrow(x) == 0) {
    warning("No rows for ", dataset_name, " ", dmr_type)
    return(NULL)
  }
  
  setorder(x, chr, start, end)
  
  ids <- make.unique(x$locus_id)
  
  vcf_dt <- data.table(
    CHROM = x$chr,
    POS = x$start,
    ID = ids,
    REF = "C",
    ALT = "T",
    QUAL = ".",
    FILTER = "PASS",
    INFO = ".",
    FORMAT = "GT"
  )
  
  for (s in samples) {
    bcol <- paste0("bin_", s)
    vcf_dt[[s]] <- bin_to_gt(x[[bcol]])
  }
  
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Synthetic genotype from binary methylation state">',
    paste(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
      paste(samples, collapse = "\t"),
      sep = "\t"
    )
  )
  
  out_vcf <- paste0(out_prefix, ".", dataset_name, ".", dmr_type, ".vcf")
  
  con <- file(out_vcf, open = "wt")
  writeLines(vcf_header, con)
  write.table(
    vcf_dt,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  close(con)
  
  if (bgzip && nzchar(Sys.which("bgzip"))) {
    system2("bgzip", args = c("-f", out_vcf))
  }
}

for (ds in c("ONT", "BIseq")) {
  for (dmr in dmr_types) {
    write_mvcf(
      df,
      dataset_name = ds,
      dmr_type = dmr,
      out_prefix = file.path(outDir, "02.4_dmr_binary_states"),
      bgzip = TRUE
    )
  }
}

# =========================================================
# 07) Raw methylation heatmaps
# =========================================================
prepare_meth_for_heatmap <- function(mat) {
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  mat <- mat[row_sd > 0 & !is.na(row_sd), , drop = FALSE]
  if (nrow(mat) == 0) return(mat)
  mat[is.na(mat)] <- -0.1
  mat
}

plot_raw_heatmap <- function(dt, dataset_name, dmr_type, heatmap_n = 5000) {
  x <- copy(dt[dataset == dataset_name & DMR == dmr_type])
  if (nrow(x) == 0) return(NULL)
  
  x[, ID := make.unique(locus_id)]
  
  mat <- as.matrix(x[, ..samples])
  rownames(mat) <- x$ID
  
  rank_dt <- data.table(
    dataset = dataset_name,
    DMR = dmr_type,
    ID = rownames(mat),
    meth_var = apply(mat, 1, var, na.rm = TRUE),
    meth_range = apply(mat, 1, safe_range),
    chr = x$chr,
    start = x$start,
    end = x$end,
    dmr_id = x$dmr_id
  )
  
  rank_dt <- rank_dt[!is.na(meth_var)]
  
  fwrite(
    rank_dt[order(-meth_var)],
    file = file.path(outDir, paste0("02.5_", dataset_name, "_", dmr_type, "_raw_methylation_ranked_by_variability.tsv")),
    sep = "\t",
    quote = FALSE
  )
  
  if (nrow(rank_dt) == 0) return(NULL)
  
  top_ids <- rank_dt[order(-meth_var)][seq_len(min(heatmap_n, .N)), ID]
  mat_top <- mat[top_ids, , drop = FALSE]
  mat_top <- prepare_meth_for_heatmap(mat_top)
  
  if (nrow(mat_top) == 0) return(NULL)
  
  ann_col <- data.frame(Group = ifelse(grepl("^G0", samples), "G0", "G5"))
  rownames(ann_col) <- samples
  p<-pheatmap(
    mat_top,
    color = colorRampPalette(c("grey90", "#4575b4", "#fee090", "#d73027"))(100),
    breaks = seq(-0.1, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    annotation_col = ann_col,
    show_rownames = FALSE,
    show_colnames = TRUE,
    main = paste0(dataset_name, " ", dmr_type, " raw methylation"),
    border_color = NA
  )
  pdf(file.path(outDir, paste0("02.5_", dataset_name, "_", dmr_type, "_raw_methylation_heatmap_top_variable.pdf")),width = 7,    height = 10  )
  print(p)
  dev.off()
}

for (ds in c("ONT", "BIseq")) {
  for (dmr in dmr_types) {
    plot_raw_heatmap(df, ds, dmr, HEATMAP_N)
  }
}

# =========================================================
# 08) Binary-state heatmaps
# =========================================================
prepare_binary_for_heatmap <- function(mat) {
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  mat <- mat[row_sd > 0 & !is.na(row_sd), , drop = FALSE]
  if (nrow(mat) == 0) return(mat)
  mat[is.na(mat)] <- -1
  mat
}

plot_binary_heatmap <- function(dt, dataset_name, dmr_type, heatmap_n = 5000) {
  x <- copy(dt[dataset == dataset_name & DMR == dmr_type])
  if (nrow(x) == 0) return(NULL)
  
  bin_cols <- paste0("bin_", samples)
  x[, ID := make.unique(locus_id)]
  
  mat <- as.matrix(x[, ..bin_cols])
  colnames(mat) <- samples
  rownames(mat) <- x$ID
  
  rank_dt <- data.table(
    dataset = dataset_name,
    DMR = dmr_type,
    ID = rownames(mat),
    binary_var = apply(mat, 1, var, na.rm = TRUE),
    binary_range = apply(mat, 1, safe_range),
    chr = x$chr,
    start = x$start,
    end = x$end,
    dmr_id = x$dmr_id
  )
  
  rank_dt <- rank_dt[!is.na(binary_var)]
  
  fwrite(
    rank_dt[order(-binary_var)],
    file = file.path(outDir, paste0("02.6_", dataset_name, "_", dmr_type, "_binary_state_ranked_by_variability.tsv")),
    sep = "\t",
    quote = FALSE
  )
  
  if (nrow(rank_dt) == 0) return(NULL)
  
  top_ids <- rank_dt[order(-binary_var)][seq_len(min(heatmap_n, .N)), ID]
  mat_top <- mat[top_ids, , drop = FALSE]
  mat_top <- prepare_binary_for_heatmap(mat_top)
  
  if (nrow(mat_top) == 0) return(NULL)
  
  ann_col <- data.frame(Group = ifelse(grepl("^G0", samples), "G0", "G5"))
  rownames(ann_col) <- samples
  p<-pheatmap(
    mat_top,
    color = c("grey90", "#2C7BB6", "#D7191C"),
    breaks = c(-1.5, -0.5, 0.5, 1.5),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    annotation_col = ann_col,
    show_rownames = FALSE,
    show_colnames = TRUE,
    main = paste0(dataset_name, " ", dmr_type, " binary states"),
    legend_breaks = c(-1, 0, 1),
    legend_labels = c("./.", "0/0", "1/1"),
    border_color = NA
  )
  
  pdf(file.path(outDir, paste0("02.6_", dataset_name, "_", dmr_type, "_binary_state_heatmap_top_variable.pdf")), width = 7,    height = 10  )
  print(p)
  dev.off()
}

for (ds in c("ONT", "BIseq")) {
  for (dmr in dmr_types) {
    plot_binary_heatmap(df, ds, dmr, HEATMAP_N)
  }
}

# =========================================================
# 09) Diagnostic plots, both datasets
# =========================================================
diag_long <- melt(df,
  id.vars = c("dataset", "chr", "start", "end", "dmr_id", "locus_id", "DMR"),
  measure.vars = samples,
  variable.name = "sample",
  value.name = "meth")

diag_long[, Group := ifelse(grepl("^G0", sample), "G0", "G5")]

# Missingness
missing_summary <- diag_long[, .(n_total = .N, n_observed = sum(!is.na(meth)),
  missing_frac = mean(is.na(meth))), by = .(dataset, DMR, sample)]

fwrite(missing_summary, file = file.path(outDir, "02.7_missingness_per_sample.tsv"), sep = "\t",  quote = FALSE)

p_missing <- ggplot(missing_summary, aes(x = sample, y = missing_frac, fill = dataset)) +
  geom_col(position = "dodge") +
  facet_grid(DMR ~ dataset, scales = "free_y") +
  theme_minimal() +
  labs(title = "Missing methylation values per sample", x = "Sample", y = "Fraction missing") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(outDir, "02.7_missingness_per_sample.pdf"),p_missing,  width = 10,  height = 6)

# Methylation density
p_density <- ggplot(diag_long[!is.na(meth)], aes(x = meth, color = sample)) +
  geom_density(linewidth = 0.8) +
  facet_grid(DMR ~ dataset, scales = "free_y") +
  theme_minimal() +
  labs(title = "Methylation level distribution per sample", x = "Methylation level",y = "Density"  )

ggsave(file.path(outDir, "02.8_methylation_density_per_sample.pdf"), p_density, width = 10,  height = 6)

# Sample deviation from row mean
df_diag <- copy(df)
meth_mat <- as.matrix(df_diag[, ..samples])
row_mean <- rowMeans(meth_mat, na.rm = TRUE)

for (s in samples) {
  df_diag[, paste0("dev_", s) := abs(get(s) - row_mean)]
}

dev_cols <- paste0("dev_", samples)

dev_long <- melt(df_diag,
  id.vars = c("dataset", "chr", "start", "end", "dmr_id", "locus_id", "DMR"),
  measure.vars = dev_cols,
  variable.name = "sample",
  value.name = "abs_dev")

dev_long[, sample := sub("^dev_", "", sample)]

dev_summary <- dev_long[, .(mean_abs_dev = mean(abs_dev, na.rm = TRUE), 
                            median_abs_dev = median(abs_dev, na.rm = TRUE)), by = .(dataset, DMR, sample)]

fwrite(dev_summary, file = file.path(outDir, "02.9_sample_deviation_from_row_mean.tsv"),sep = "\t",  quote = FALSE)

p_dev <- ggplot(dev_summary, aes(x = sample, y = mean_abs_dev, fill = dataset)) +
  geom_col(position = "dodge") +
  facet_grid(DMR ~ dataset, scales = "free_y") +
  theme_minimal() +
  labs(title = "Mean absolute methylation deviation from DMR row mean",
    x = "Sample",    y = "Mean absolute deviation"  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(outDir, "02.9_sample_deviation_from_row_mean.pdf"),p_dev,  width = 10,  height = 6)

# G0-C6 divergent loci
df_diag[, row_mean_no_G0C6 := rowMeans(.SD, na.rm = TRUE), .SDcols = setdiff(samples, "G0-C6")]
df_diag[, G0C6_effect := abs(`G0-C6` - row_mean_no_G0C6)]

g0c6_top <- df_diag[!is.na(G0C6_effect)]
g0c6_top <- g0c6_top[
  order(dataset, DMR, -G0C6_effect),
  .SD[seq_len(min(1000, .N))],
  by = .(dataset, DMR)]

fwrite(g0c6_top, file = file.path(outDir, "02.10_top_G0C6_divergent_DMRs.tsv"), sep = "\t",  quote = FALSE)

p_g0c6 <- ggplot(g0c6_top, aes(x = start, fill = dataset)) +
  geom_histogram(bins = 100) +
  facet_grid(dataset + DMR ~ chr, scales = "free_x") +
  theme_minimal() +
  labs(title = "Genomic distribution of top G0-C6 divergent DMRs", x = "Position",y = "Count")

ggsave(file.path(outDir, "02.10_G0C6_divergent_DMRs_genomic_distribution.pdf"), p_g0c6,  width = 14,  height = 8)