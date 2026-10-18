suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# Paths
# =========================

workdir <- "/mnt/ssd123/vibanez/19_ont-mapping/ab_run"
contig_tab <- file.path(workdir, "polished_contigs_used_in_chr1to12.tsv")
sl5_fai  <- "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SL5.0.fasta.fai"
sl25_fai <- "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta.fai"
outdir <- "/mnt/disk2/vibanez/otherAnalysis/review_denovo-assembly/contig_size_plots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

MIN_CONTIG_BP <- 1e6

# =========================
# Cervil de novo chromosome sizes
# from polished contigs used by RagTag
# =========================

D <- fread(contig_tab)

D[, chr := gsub("_RagTag", "", ragtag_chr)]
D <- D[chr %in% as.character(1:12)]
D[, chr := factor(chr, levels = as.character(1:12))]

cervil <- D[, .(
  size_bp = sum(placed_bp),
  n_contigs = .N,
  n_large_contigs = sum(placed_bp >= MIN_CONTIG_BP)
), by = chr]

cervil[, genome := "Cervil de novo"]
cervil[, size_Mb := size_bp / 1e6]

# =========================
# SL5 chromosome sizes
# =========================

sl5 <- fread(sl5_fai, header = FALSE)
setnames(sl5, c("seqid", "length", "offset", "linebases", "linewidth"))

sl5 <- sl5[seqid %in% as.character(1:12)]
sl5[, chr := factor(seqid, levels = as.character(1:12))]
sl5[, genome := "SL5"]
sl5[, size_bp := length]
sl5[, size_Mb := size_bp / 1e6]
sl5[, `:=`(n_contigs = NA_integer_, n_large_contigs = NA_integer_)]

# =========================
# SL2.5 chromosome sizes
# =========================

sl25 <- fread(sl25_fai, header = FALSE)
setnames(sl25, c("seqid", "length", "offset", "linebases", "linewidth"))

sl25[, chr_num := sub(".*ch0*([0-9]+)$", "\\1", seqid)]
sl25 <- sl25[chr_num %in% as.character(1:12)]
sl25[, chr := factor(chr_num, levels = as.character(1:12))]
sl25[, genome := "SL2.5"]
sl25[, size_bp := length]
sl25[, size_Mb := size_bp / 1e6]
sl25[, `:=`(n_contigs = NA_integer_, n_large_contigs = NA_integer_)]

# =========================
# Combine
# =========================

plot_dt <- rbindlist(list(
  sl25[, .(chr, genome, size_bp, size_Mb, n_contigs, n_large_contigs)],
  sl5[,  .(chr, genome, size_bp, size_Mb, n_contigs, n_large_contigs)],
  cervil[, .(chr, genome, size_bp, size_Mb, n_contigs, n_large_contigs)]
), use.names = TRUE)

plot_dt[, genome := factor(
  genome,
  levels = c("SL2.5", "SL5", "Cervil de novo")
)]

# Save table
fwrite(
  plot_dt,
  file = file.path(outdir, "SL25_SL5_Cervil_denovo_chr_sizes.tsv"),
  sep = "\t"
)

# =========================
# Plot
# =========================
pd <- position_dodge(width = 0.8)
p <- ggplot(plot_dt, aes(x = chr, y = size_Mb, fill = genome)) +
  geom_col( position = pd, width = 0.72, color = "grey30", linewidth = 0.2) +
  geom_text(data = plot_dt[genome == "Cervil de novo"],
    aes(label = paste0("n=", n_large_contigs)),
    position = pd, vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("SL2.5" = "#666666", "SL5" = "#BDBDBD","Cervil de novo" = "#1B9E77")) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Chromosome", y = "Size (Mb)", fill = "Genome",
    title = "Chromosome sizes in SL2.5, SL5 and Cervil de novo assembly",
    subtitle = "n = number of contigs") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right" )

ggsave(filename = file.path(outdir, "SL25_SL5_Cervil_denovo_chr_sizes_grouped_bars.pdf"),
plot = p,  width = 10,  height = 5)
