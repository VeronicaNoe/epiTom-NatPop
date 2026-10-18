suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(gridExtra)
})


RESULTS_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/results"
PLOTS_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/plots"
#mkdir -p "${COMP_DIR}" "${RESULTS_DIR}" "${PLOTS_DIR}"

WINDOW_BP="${WINDOW_BP:-1000}"
SNP_VCF="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants/cervil_SNPs.GQ20.DP10.vcf.gz"
SV_VCF="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants/cervil_SV.PASS.SUP5.vcf.gz"

root <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants"
sample_cols <- c("G0-C3","G0-C6","G0-C9","G5-C1","G5-C2","G5-C3")
snp_eigenvec_file <- file.path(root, "cervil-SNPs-ONT_pca.eigenvec")
snp_eigenval_file <- file.path(root, "cervil-SNPs-ONT_pca.eigenval")
snp_mdist_file    <- file.path(root, "cervil-SNPs-ONT_dendro.mdist.gz")
snp_mdist_id_file <- file.path(root, "cervil-SNPs-ONT_dendro.mdist.id")

sv_eigenvec_file  <- file.path(root, "cervil-SV-ONT_pca.eigenvec")
sv_eigenval_file  <- file.path(root, "cervil-SV-ONT_pca.eigenval")
sv_mdist_file     <- file.path(root, "cervil-SV-ONT_dendro.mdist.gz")
sv_mdist_id_file  <- file.path(root, "cervil-SV-ONT_dendro.mdist.id")


snp_vcf <- fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'",
  shQuote(SNP_VCF)))

sv_vcf <- fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\t%INFO/SVTYPE\\t%INFO/SVLEN[\\t%GT]\\n'",
  shQuote(SV_VCF)))

samples <- c("G0-C3","G0-C6","G0-C9","G5-C1","G5-C2","G5-C3")

setnames(snp_vcf, c("CHROM","POS","REF","ALT", samples))
setnames(sv_vcf,  c("CHROM","POS","SVTYPE","SVLEN", samples))

snp_vcf[, variant_type := ifelse(nchar(REF)==1 & nchar(ALT)==1, "SNP", "INDEL")]
gt_mat <- as.matrix(snp_vcf[, ..samples])

all_missing <- rowSums(
  is.na(gt_mat) | gt_mat %in% c("./.", ".|.")
) == length(samples)

snp_vcf_nonmissing <- snp_vcf[!all_missing]

# -----------------------------
# Variant counts per sample/type
# -----------------------------
to_binary <- function(x) {
  ifelse(is.na(x) | x %in% c("./.", ".|."), NA,
         ifelse(x %in% c("0/0", "0|0"), 0, 1))
}

# SNP / small indel table
snp_long <- melt(snp_vcf_nonmissing, id.vars = c("CHROM", "POS", "REF", "ALT", "variant_type"),
                 measure.vars = samples,  variable.name = "sample",   value.name = "GT")

snp_long[, var_class := fifelse(variant_type == "SNP", "SNP", "small_indel")]
snp_long[, present := to_binary(GT)]

# SV table
sv_vcf[, SVLEN := as.numeric(SVLEN)]

sv_long <- melt(sv_vcf, id.vars = c("CHROM", "POS", "SVTYPE", "SVLEN"), measure.vars = samples,
  variable.name = "sample",  value.name = "GT")

sv_long[, var_class := SVTYPE]
sv_long[, present := to_binary(GT)]

# combine
var_long <- rbindlist(list(
  snp_long[, .(sample, var_class, present)],
  sv_long[,  .(sample, var_class, present)]), fill = TRUE)

# count variants present per sample
count_dt <- var_long[present == 1, .N, by = .(sample, var_class)]

# make all combinations explicit
count_dt <- CJ(sample = samples, var_class = unique(var_long$var_class))[
  count_dt, on = .(sample, var_class)
]
count_dt[is.na(N), N := 0]

count_dt[, total := sum(N), by = sample]
count_dt[, frac := N / total]

count_dt
p <- ggplot(count_dt, aes(x = sample, y = frac, fill = var_class)) +
  geom_col(width = 0.75) +
  geom_text( aes(label = ifelse(N > 0, N, "")), position = position_stack(vjust = 0.5),size = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Variant class proportion", fill = "Variant type", title = "Variant composition per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PLOTS_DIR, "01.01_variant_composition_per_sample_100pct.pdf"),p,width = 7,height = 5)

fwrite(count_dt, file.path(RESULTS_DIR, "01.01_variant_counts_per_sample.tsv"), sep = "\t")
##
# -----------------------------
# Heatmap: samples x variants using real genotypes
# -----------------------------
# Add variant IDs before combining
snp_long[, variant_id := paste(CHROM, POS, REF, ALT, sep = ":")]
sv_long[,  variant_id := paste(CHROM, POS, SVTYPE, SVLEN, sep = ":")]

heat_dt <- rbindlist(list(
  snp_long[, .(variant_id, sample, var_class, GT)],
  sv_long[,  .(variant_id, sample, var_class, GT)]
), fill = TRUE)

# Simplify genotype classes
heat_dt[, GT_class := fifelse(
  is.na(GT) | GT %in% c("./.", ".|.", "0/.", "./0", "1/.", "./1"),
  "missing",
  fifelse(GT %in% c("0/0", "0|0"),
          "0/0",
          fifelse(GT %in% c("0/1", "1/0", "0|1", "1|0"),
                  "0/1",
                  fifelse(GT %in% c("1/1", "1|1"),
                          "1/1",
                          "other")))
)]

# order variants by class and genotype pattern
pattern_dt <- dcast(
  heat_dt,
  variant_id + var_class ~ sample,
  value.var = "GT_class"
)

pattern_dt[, pattern := apply(.SD, 1, paste, collapse = ""), .SDcols = samples]

setorder(pattern_dt, var_class, pattern, variant_id)

heat_dt[, variant_id := factor(variant_id, levels = rev(pattern_dt$variant_id))]
heat_dt[, sample := factor(sample, levels = samples)]

heat_dt[, GT_class := factor(
  GT_class,
  levels = c("0/0", "0/1", "1/1", "missing", "other")
)]

p_heat <- ggplot(heat_dt, aes(x = sample, y = variant_id, fill = GT_class)) +
  geom_tile(color = NA) +
  facet_grid(var_class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("0/0" = "#2C7BB6", "0/1" = "#7B3294","1/0" = "#7B3294","1/1" = "#D7191C",
                               "missing" = "white"  ),  name = "Genotype",  drop = FALSE  ) +
  labs(x = NULL,y = NULL, title = "Genotype heatmap") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0),panel.grid = element_blank()  )

ggsave(file.path(PLOTS_DIR, "01.02_variant_genotype_heatmap.pdf"),  p_heat,  width = 6,height = 10)
#

# -----------------------------
# Clustered heatmap with pheatmap
# ----------------------------
gt_num_dt <- copy(heat_dt)

gt_num_dt[, GT_num := fifelse(
  GT_class == "0/0", 0,
  fifelse(GT_class == "0/1", 1,
          fifelse(GT_class == "1/1", 2, NA_real_))
)]

gt_wide <- dcast(
  gt_num_dt,
  variant_id + var_class ~ sample,
  value.var = "GT_num"
)

mat <- as.matrix(gt_wide[, ..samples])
rownames(mat) <- gt_wide$variant_id

# remove all-NA rows if any
mat <- mat[rowSums(!is.na(mat)) > 0, ]

# impute NA only for clustering/plotting
mat_plot <- mat
row_means <- rowMeans(mat_plot, na.rm = TRUE)
row_means[is.na(row_means)] <- 0

for (i in seq_len(nrow(mat_plot))) {
  mat_plot[i, is.na(mat_plot[i, ])] <- row_means[i]
}

# row annotation: variant class
ann_row <- data.frame(
  Variant_class = gt_wide$var_class[match(rownames(mat_plot), gt_wide$variant_id)]
)
rownames(ann_row) <- rownames(mat_plot)

# column annotation: generation/group
ann_col <- data.frame(
  Group = ifelse(grepl("^G0", samples), "G0", "G5")
)
rownames(ann_col) <- samples

pdf(file.path(PLOTS_DIR, "01.03_variant_genotype_pheatmap_clustered.pdf"),
    width = 6, height = 10)

pheatmap(
  mat_plot,
  color = c("#2C7BB6", "#7B3294", "#D7191C"),
  breaks = c(-0.5, 0.5, 1.5, 2.5),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = ann_row,
  annotation_col = ann_col,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Clustered genotype heatmap",
  legend_breaks = c(0, 1, 2),
  legend_labels = c("0/0", "0/1", "1/1"),
  border_color = NA
)

dev.off()
###
# -----------------------------
# SV length distribution
# -----------------------------
sv_len_dt <- copy(sv_vcf)
sv_len_dt[, abs_SVLEN := abs(as.numeric(SVLEN))]
sv_len_dt <- sv_len_dt[!is.na(abs_SVLEN) & abs_SVLEN > 0]

len_levels <- c("1-50", "51-100", "101-500", "501-1kb", "1-5kb",
                "5-10kb", "10-50kb", ">50kb")

sv_len_dt[, len_bin := cut(abs_SVLEN,
  breaks = c(0, 50, 100, 500, 1000, 5000, 10000, 50000, Inf),
  labels = len_levels,
  right = TRUE)]

sv_len_dt[, len_bin := factor(len_bin, levels = len_levels)]

# Count SVs by type and length bin
sv_len_counts <- sv_len_dt[, .N, by = .(SVTYPE, len_bin)]

# Make empty combinations explicit
sv_len_counts <- CJ(
  SVTYPE = unique(sv_len_dt$SVTYPE),
  len_bin = len_levels)[
  sv_len_counts,
  on = .(SVTYPE, len_bin)
]

sv_len_counts[is.na(N), N := 0]
sv_len_counts[, len_bin := factor(len_bin, levels = len_levels)]

fwrite(sv_len_counts,file.path(RESULTS_DIR, "01.04_sv_length_bins_by_type.tsv"),sep = "\t")

# -----------------------------
# Count histogram
# -----------------------------
p_sv_count <- ggplot(sv_len_counts, aes(x = len_bin, y = N, fill = SVTYPE)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(N > 0, N, "")),position = position_stack(vjust = 0.5),size = 3) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "SV length bin (bp)",y = "Number of SVs",fill = "SV type",title = "SV length distribution") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PLOTS_DIR, "01.04_sv_length_distribution_counts.pdf"),p_sv_count,width = 7,height = 5)

# -----------------------------
# Frequency histogram
# frequency = proportion of each length bin
# -----------------------------
sv_len_counts[, freq := {
  total <- sum(N)
  if (total > 0) N / total else rep(0, .N)
}, by = len_bin]
p_sv_freq <- ggplot(sv_len_counts, aes(x = len_bin, y = freq, fill = SVTYPE)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(N > 0, N, "")),position = position_stack(vjust = 0.5),size = 3) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "SV length bin (bp)",y = "Frequency within length bin",
    fill = "SV type",title = "SV length distribution frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PLOTS_DIR, "01.05_sv_length_distribution_frequency.pdf"),p_sv_freq,width = 7,height = 5)
# -----------------------------
# SV length per sample and SV type
# -----------------------------

sv_len_sample <- copy(sv_long)

sv_len_sample[, SVLEN := abs(as.numeric(SVLEN))]
sv_len_sample[, present := to_binary(GT)]

sv_len_sample <- sv_len_sample[present == 1 & !is.na(SVLEN)]

len_levels <- c("1-50", "51-100", "101-500", "501-1kb", "1-5kb","5-10kb", "10-50kb", ">50kb")

sv_len_sample[, len_bin := cut(SVLEN, 
                               breaks = c(0, 50, 100, 500, 1000, 5000, 10000, 50000, Inf),
                               labels = len_levels,right = TRUE)]

sv_len_sample[, len_bin := factor(len_bin, levels = len_levels)]

# count SVs by sample, type, and length bin
sv_len_sample_count <- sv_len_sample[, .N, by = .(sample, SVTYPE, len_bin)]

# make empty combinations explicit
sv_len_sample_count <- CJ(sample = samples,
                          SVTYPE = unique(sv_len_sample$SVTYPE),
                          len_bin = len_levels)[
                            sv_len_sample_count,
                            on = .(sample, SVTYPE, len_bin)]

sv_len_sample_count[is.na(N), N := 0]
sv_len_sample_count[, len_bin := factor(len_bin, levels = len_levels)]

p_sv_sample_len <- ggplot(sv_len_sample_count,aes(x = len_bin, y = N, fill = SVTYPE)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = ifelse(N > 0, N, "")),position = position_stack(vjust = 0.5),size = 2.8) +
  facet_wrap(~ sample, nrow = 2) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "SV length bin (bp)",y = "Number of carried SVs",
    fill = "SV type",title = "SV length distribution per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.minor = element_blank())

ggsave(file.path(PLOTS_DIR, "01.06_sv_length_distribution_per_sample.pdf"), p_sv_sample_len, width = 10,height = 6)

fwrite(sv_len_sample_count,file.path(RESULTS_DIR, "01.06_sv_length_distribution_per_sample.tsv"),sep = "\t")
# -----------------------------
# Pairwise genetic distance
# -----------------------------
# wide binary matrix: rows = variants, cols = samples
heat_bin <- copy(heat_dt)
heat_bin[, present := fifelse(GT_class == "0/0", 0,
                              fifelse(GT_class %in% c("0/1", "1/1"), 1, NA_real_))]

bin_wide <- dcast(heat_bin, variant_id + var_class ~ sample, value.var = "present")

pairwise_dt <- rbindlist(lapply(unique(bin_wide$var_class), function(vc) {
  tmp <- bin_wide[var_class == vc]
  
  rbindlist(lapply(samples, function(s1) {
    rbindlist(lapply(samples, function(s2) {
      x <- tmp[[s1]]
      y <- tmp[[s2]]
      valid <- !is.na(x) & !is.na(y)
      data.table(
        var_class = vc,
        sample1 = s1,
        sample2 = s2,
        n_compared = sum(valid),
        n_different = sum(x[valid] != y[valid]),
        frac_different = sum(x[valid] != y[valid]) / sum(valid)
      )
    }))
  }))
}))

pairwise_dt[, label := ifelse(
  is.na(frac_different),
  "NA",
  sprintf("%.2f", frac_different)
)]

p_dist <- ggplot(pairwise_dt, aes(sample1, sample2, fill = frac_different)) +
  geom_tile(color = "grey90", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3) +
  facet_wrap(~ var_class) +
  scale_fill_gradient(
    low = "white",
    high = "#D7191C",
    na.value = "grey80",
    name = "Fraction\ndifferent",
    limits = c(0, 1)) +
  labs(x = NULL,y = NULL,title = "Pairwise genetic distance by variant type") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "01.07_pairwise_distance_by_variant_type.pdf"),p_dist,width = 9,height = 6)
fwrite(pairwise_dt, file.path(RESULTS_DIR, "01.06_pairwise_distance_by_variant_type.tsv"), sep = "\t")
##

# -----------------------------
# PCA plotting function
# -----------------------------
plot_plink_pca <- function(eigenvec_file, eigenval_file, prefix, outdir) {
  
  ev <- fread(eigenvec_file, header = FALSE)
  colnames(ev) <- c("FID", "IID", paste0("PC", seq_len(ncol(ev) - 2)))
  
  eig <- scan(eigenval_file, quiet = TRUE)
  var_exp <- 100 * eig / sum(eig)
  
  ev[, group := ifelse(grepl("^G0", IID), "G0", "G5")]
  ev[, IID := factor(IID, levels = sample_cols)]
  
  p <- ggplot(ev, aes(PC1, PC2, label = IID, color = group)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    theme_bw() +
    labs(
      x = sprintf("PC1 (%.2f%%)", var_exp[1]),
      y = sprintf("PC2 (%.2f%%)", var_exp[2]),
      title = paste(prefix, "PCA"),
      color = "Group"
    )
  
  fwrite(ev, file.path(outdir, paste0("01.08_",prefix, "_PCA_scores.tsv")), sep = "\t")
  
  ggsave(file.path(outdir, paste0("01.08_",prefix, "_PCA.pdf")),p,width = 6,    height = 5  )
  
  return(p)
}

# -----------------------------
# Dendrogram plotting function
# -----------------------------
plot_plink_tree <- function(mdist_file, mdist_id_file, prefix, outdir) {
  
  ids <- fread(mdist_id_file, header = FALSE)
  md  <- as.matrix(fread(mdist_file, header = FALSE))
  
  samples <- ids$V2
  rownames(md) <- samples
  colnames(md) <- samples
  
  fwrite(
    as.data.table(md, keep.rownames = "sample"),
    file.path(outdir, paste0(prefix, "_distance_matrix.tsv")),
    sep = "\t"
  )
  
  hc <- hclust(as.dist(md), method = "average")
  
  pdf(file.path(outdir, paste0("01.09_",prefix, "_dendrogram.pdf")), width = 7, height = 5)
  plot(hc, main = paste("01.09_",prefix, "distance tree"), xlab = "", sub = "")
  dev.off()
  
  return(hc)
}

# -----------------------------
# Run SNP and SV PCA/dendrogram
# -----------------------------
p_snp <- plot_plink_pca(snp_eigenvec_file, snp_eigenval_file, "SNP",  PLOTS_DIR)

p_sv <- plot_plink_pca(sv_eigenvec_file,sv_eigenval_file,  "SV",  PLOTS_DIR)

hc_snp <- plot_plink_tree(snp_mdist_file,  snp_mdist_id_file,  "SNP",  PLOTS_DIR)

hc_sv <- plot_plink_tree(sv_mdist_file, sv_mdist_id_file,"SV",  PLOTS_DIR)

# -----------------------------
# 4) Combined plots
# -----------------------------

pdf(file.path(outdir, "SNP_SV_PCA_combined.pdf"), width = 12, height = 5)
grid.arrange(p_snp, p_sv, ncol = 2)
dev.off()
pdf(file.path(outdir, "SNP_SV_dendrogram_combined.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))
plot(hc_snp, main = "SNP distance tree", xlab = "", sub = "")
plot(hc_sv,  main = "SV distance tree",  xlab = "", sub = "")
dev.off()
par(mfrow = c(1, 1))
