suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
root <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants"

FAI <- "/mnt/ssd123/vibanez/19_ont-mapping/ab_run/cervil.polished.fasta.fai"

SNP_VCF <- file.path(root, "cervil_SNPs.GQ20.DP10.vcf.gz")
SV_VCF  <- file.path(root, "cervil_SV.PASS.SUP5.vcf.gz")

RESULTS_DIR <- file.path(root, "results")
PLOTS_DIR <- file.path(root, "plots")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

samples <- c("G0-C3","G0-C6","G0-C9","G5-C1","G5-C2","G5-C3")

# -----------------------------
# Parameters
# -----------------------------
TOP_N <- 12
MIN_CONTIG_BP <- 10000000   # adjust after checking contig ranking

# "segregating" = different genotypes among samples
# "alt_present" = any sample carries ALT
MODE <- "segregating"

# -----------------------------
# 1. Contig lengths
# -----------------------------
fai <- fread(FAI, header = FALSE)
fai <- fai[, .(CHROM = V1, length = V2)]

setorder(fai, -length)
fai[, rank := .I]

fwrite(fai, file.path(RESULTS_DIR, "02.01_contig_lengths_ranked.tsv"),sep = "\t")

# fai already has: CHROM, length, rank
fai[, length_Mb := length / 1e6]

fai_top20 <- fai[order(-length)][1:20]

fai_top20[, CHROM := factor(CHROM, levels = rev(CHROM))]

p_top20 <- ggplot(fai_top20, aes(x = CHROM, y = length_Mb)) +
  geom_col(fill = "grey60") +
  geom_text(aes(label = sprintf("%.1f Mb", length_Mb)), hjust = -0.1,size = 3  ) +
  coord_flip() +
  expand_limits(y = max(fai_top20$length_Mb) * 1.12) +
  labs( x = NULL, y = "Contig length (Mb)", title = "Top 20 longest contigs") +
  theme_bw()
ggsave(filename = file.path(PLOTS_DIR, "04.01_contig_length_rank_plot.pdf"), plot = p_top20 ,  width = 7,  height = 5)

long_contigs <- fai[rank <= TOP_N & length >= MIN_CONTIG_BP]
fwrite(long_contigs, file.path(RESULTS_DIR, "02.02_chromosome_like_contigs.tsv"), sep = "\t")

long_contig_bed <- long_contigs[, .(CHROM, start = 0,end = length)]

fwrite(long_contig_bed,file.path(RESULTS_DIR, "02.02_chromosome_like_contigs.bed"),
  sep = "\t",  col.names = FALSE)

keep_contigs <- long_contigs$CHROM

cat("Keeping contigs:\n")
print(long_contigs)

# -----------------------------
# Helper functions
# -----------------------------
norm_gt_vec <- function(x) {
  x <- gsub("\\|", "/", x)
  x[is.na(x) | grepl("\\.", x)] <- NA_character_
  
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  
  out[ok] <- vapply(
    strsplit(x[ok], "/", fixed = TRUE),
    function(a) paste(sort(a), collapse = "/"),
    character(1)
  )
  
  out
}

keep_variant_rows <- function(dt, samples, mode = c("segregating", "alt_present")) {
  mode <- match.arg(mode)
  
  gt_norm <- as.data.table(lapply(dt[, ..samples], norm_gt_vec))
  
  if (mode == "segregating") {
    keep <- apply(gt_norm, 1, function(x) {
      x <- x[!is.na(x)]
      length(unique(x)) >= 2
    })
  } else {
    keep <- apply(gt_norm, 1, function(x) {
      x <- x[!is.na(x)]
      any(grepl("[1-9]", x))
    })
  }
  
  keep
}

max_alt_len <- function(alt) {
  vapply(strsplit(alt, ","), function(x) max(nchar(x)), numeric(1))
}

parse_num1 <- function(x) {
  suppressWarnings(as.numeric(sub(",.*", "", x)))
}

write_bed <- function(dt, file) {
  fwrite(dt, file, sep = "\t", col.names = FALSE)
}

merge_bed_intervals <- function(bed) {
  if (nrow(bed) == 0) {
    return(bed[, .(CHROM, start, end)])
  }
  
  x <- copy(bed[, .(CHROM, start, end)])
  setorder(x, CHROM, start, end)
  
  x[, prev_max_end := shift(cummax(end), fill = -Inf), by = CHROM]
  x[, grp := cumsum(start > prev_max_end), by = CHROM]
  
  out <- x[, .(
    start = min(start),
    end = max(end),
    n_intervals = .N
  ), by = .(CHROM, grp)]
  
  out[, grp := NULL]
  setorder(out, CHROM, start, end)
  out
}

# -----------------------------
# 2. SNP / small INDEL variants
# -----------------------------
snp <- fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'",shQuote(SNP_VCF)))

setnames(snp, c("CHROM", "POS", "REF", "ALT", samples))

snp <- snp[CHROM %in% keep_contigs]

snp[, is_variable := keep_variant_rows(.SD, samples, MODE), .SDcols = names(snp)]

snp_var <- snp[is_variable == TRUE]

snp_var[, variant_type := ifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNP", "small_indel")]

snp_var[, start := POS - 1]
snp_var[, end := start + pmax(nchar(REF), max_alt_len(ALT), 1)]
snp_var[, name := paste(CHROM, POS, REF, ALT, sep = ":")]

snp_bed <- snp_var[, .(CHROM, start,  end,  name,  score = 0,  strand = ".",  variant_type)]

write_bed( snp_bed, file.path(RESULTS_DIR, "02.03_variable_SNPs_INDELs.bed"))

# -----------------------------
# 3. SV variants
# -----------------------------
sv <- fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN[\\t%GT]\\n'",shQuote(SV_VCF)))

setnames(sv, c("CHROM", "POS", "END", "SVTYPE", "SVLEN", samples))

sv <- sv[CHROM %in% keep_contigs]

sv[, is_variable := keep_variant_rows(.SD, samples, MODE), .SDcols = names(sv)]

sv_var <- sv[is_variable == TRUE]

sv_var[, END_num := parse_num1(END)]
sv_var[, SVLEN_num := abs(parse_num1(SVLEN))]

sv_var[, start := POS - 1]

sv_var[, end := fifelse(!is.na(END_num) & END_num > POS,
                        END_num,fifelse(!is.na(SVLEN_num) & SVLEN_num > 0,
                                        POS + SVLEN_num - 1, POS))]

sv_var[end <= start, end := start + 1]

sv_var[, name := paste(CHROM, POS, SVTYPE, SVLEN, sep = ":")]

sv_bed <- sv_var[, .(CHROM,  start,  end,  name,  score = 0,strand = ".",  variant_type = SVTYPE)]

write_bed(sv_bed,file.path(RESULTS_DIR, "02.04_variable_SVs.bed"))

# -----------------------------
# 4. Combined variant mask
# -----------------------------
all_bed <- rbindlist(list(snp_bed, sv_bed), fill = TRUE)
setorder(all_bed, CHROM, start, end)

write_bed(  all_bed,  file.path(RESULTS_DIR, "02.05_variable_variants.all.bed"))

all_bed_merged <- merge_bed_intervals(all_bed)

write_bed(all_bed_merged, file.path(RESULTS_DIR, "02.06_variable_variants.all.merged.bed"))

# -----------------------------
# 4b. Combined variant genotype table
# -----------------------------

snp_gt <- snp_var[, c("CHROM", "start", "end", "name", "variant_type","POS", "REF", "ALT",  samples), with = FALSE]

setnames(snp_gt, "name", "variant_id")

snp_gt[, `:=`(variant_class = "SNP_or_small_indel", END = NA_character_,   SVTYPE = NA_character_,  SVLEN = NA_character_)]

sv_gt <- sv_var[, c("CHROM", "start", "end", "name",  "SVTYPE", "POS", "END", "SVLEN",  samples), with = FALSE]

setnames(sv_gt, "name", "variant_id")

sv_gt[, variant_class := "SV"]
sv_gt[, variant_type := SVTYPE]
sv_gt[, REF := NA_character_]
sv_gt[, ALT := NA_character_]

gt_cols <- c("CHROM", "start", "end",  "variant_id", "variant_class", "variant_type",  "POS", "REF", "ALT", "END", "SVTYPE", "SVLEN",  samples)

all_gt <- rbindlist(list(snp_gt[, ..gt_cols],    sv_gt[, ..gt_cols]  ),  fill = TRUE)

setorder(all_gt, CHROM, start, end, variant_id)

# normalized genotype columns
norm_cols <- paste0(samples, "_normGT")
all_gt[, (norm_cols) := lapply(.SD, norm_gt_vec), .SDcols = samples]

fwrite(all_gt, file.path(RESULTS_DIR, "02.05_variable_variants.all.genotypes.tsv"), sep = "\t")
# -----------------------------
# 5. Summary
# -----------------------------
summary_dt <- rbindlist(list(
  data.table(type = "SNP_or_small_indel", n = nrow(snp_bed)),
  data.table(type = "SV", n = nrow(sv_bed)),
  data.table(type = "all_exact_intervals", n = nrow(all_bed))))

fwrite(summary_dt, file.path(RESULTS_DIR, "02.07_variable_variant_mask_summary.tsv"),sep = "\t")

print(summary_dt)
