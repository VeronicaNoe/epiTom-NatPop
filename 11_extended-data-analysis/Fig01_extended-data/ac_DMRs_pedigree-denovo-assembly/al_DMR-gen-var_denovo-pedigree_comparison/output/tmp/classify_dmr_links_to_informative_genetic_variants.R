#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly"

METH_DIR <- file.path(BASE, "ai_DMR-classification/bd_merged-methylation")
GEN_DIR  <- file.path(BASE, "ak_genetic-variants")
COMP_DIR <- file.path(BASE, "al_DMR-gen-var_denovo-pedigree_comparison")
RESULTS_DIR <- file.path(BASE, "results")
PLOTS_DIR   <- file.path(BASE, "plots")

SNP_VCF <- file.path(GEN_DIR, "cervil_SNPs.GQ20.DP10.vcf.gz")
SV_VCF  <- file.path(GEN_DIR, "cervil_SV.PASS.SUP5.vcf.gz")


WINDOW_BP <- 1000L
MIN_SUPPORT <- 1
MIN_MATCH_FRAC <- 0.15

run_query <- function(cmd) {
  fread(cmd = cmd, sep = "\t", header = FALSE, data.table = TRUE, fill = TRUE)
}

gt_to_present <- function(x) {
  x <- as.character(x)
  fifelse(
    is.na(x) | x %in% c(".", "./.", ".|.", "0/.", "./0", "1/.", "./1"),
    NA_integer_,
    fifelse(x %in% c("0/0", "0|0"), 0L, 1L)
  )
}

gt_to_meth_state <- function(x) {
  x <- as.character(x)
  fifelse(
    is.na(x) | x %in% c(".", "./.", ".|.", "0/.", "./0", "1/.", "./1"),
    NA_integer_,
    fifelse(x %in% c("0/0", "0|0"), 0L,
            fifelse(x %in% c("1/1", "1|1"), 1L, NA_integer_))
  )
}

is_informative_pattern <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 2) return(FALSE)
  length(unique(v)) >= 2
}

calc_pair_support <- function(meth, var) {
  ok <- !is.na(meth) & !is.na(var)
  n <- sum(ok)

  if (n == 0) {
    return(list(
      n_compared = 0L,
      direct_support = 0L,
      inverse_support = 0L,
      best_support = 0L,
      best_match_frac = NA_real_,
      direction = NA_character_
    ))
  }

  direct_support <- sum(meth[ok] == var[ok])
  inverse_support <- sum(meth[ok] != var[ok])

  if (direct_support >= inverse_support) {
    best_support <- direct_support
    direction <- "direct"
  } else {
    best_support <- inverse_support
    direction <- "inverse"
  }

  list(
    n_compared = n,
    direct_support = direct_support,
    inverse_support = inverse_support,
    best_support = best_support,
    best_match_frac = best_support / n,
    direction = direction
  )
}

load_dmr_one <- function(dmr_type) {
  message("Loading DMR states: ", dmr_type)

  dmr_vcf <- file.path(METH_DIR, paste0("00.1_dmr_binary_states.", dmr_type, ".vcf.gz"))
  meth_bed <- file.path(METH_DIR, paste0(dmr_type, "_all.merged.methylation.bed"))
  if (!file.exists(dmr_vcf)) stop("Missing DMR VCF: ", dmr_vcf)
  if (!file.exists(meth_bed)) stop("Missing methylation BED: ", meth_bed)

  q <- paste("bcftools query -f '%CHROM\\t%POS0\\t%ID[\\t%GT]\\n'",shQuote(dmr_vcf)  )

  st <- run_query(q)
  setnames(st, c("chr", "start", "DMRid", samples))
  st[, start := as.integer(start)]

  for (s in samples) {
    st[, paste0("meth_", s) := gt_to_meth_state(get(s))]
  }
  st <- st[, c("chr", "start", "DMRid", paste0("meth_", samples)), with = FALSE]
  bed <- fread(meth_bed, header = FALSE, sep = "\t", data.table = TRUE, fill = TRUE)
  setnames(bed, 1:3, c("chr", "start", "end"))
  bed <- unique(bed[, .(chr, start = as.integer(start), end = as.integer(end))])

  x <- merge(st, bed, by = c("chr", "start"), all.x = TRUE)
  x[is.na(end), end := start + 1L]

  x[, DMR := dmr_type]
  setcolorder(x, c("chr", "start", "end", "DMRid", "DMR", paste0("meth_", samples)))
  x
}

samples <- c("G0-C3", "G0-C6", "G0-C9", "G5-C1", "G5-C2", "G5-C3")
dmr_types <- c("C-DMR", "CG-DMR")

# -------------------------------------------------------
# load DMRs
# -------------------------------------------------------
dmr_dt <- rbindlist(lapply(dmr_types, load_dmr_one), use.names = TRUE, fill = TRUE)
dmr_dt <- unique(dmr_dt)

fwrite(dmr_dt, file.path(COMP_DIR, "01_dmr_binary_states_full_intervals.tsv"), sep = "\t")

# -------------------------------------------------------
# load SNPs
# -------------------------------------------------------
message("Loading SNP/small-indel VCF...")

snp_q <- paste("bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n'",shQuote(SNP_VCF))

snp <- run_query(snp_q)
setnames(snp, c("chr", "pos", "variant_id", "REF", "ALT", samples))
snp[, pos := as.integer(pos)]
snp[variant_id == "." | is.na(variant_id), variant_id := paste(chr, pos, REF, ALT, sep = ":")]

snp[, var_class := fifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNP", "small_indel")]
snp[, start := pos - 1L]
snp[, end := start + pmax(nchar(REF), nchar(ALT), 1L)]

for (s in samples) {
  snp[, paste0("var_", s) := gt_to_present(get(s))]
}

snp[, informative := apply(.SD, 1, is_informative_pattern), .SDcols = paste0("var_", samples)]
snp_inf <- snp[informative == TRUE]

snp_inf <- snp_inf[, c(
  "chr", "start", "end", "variant_id", "var_class", "REF", "ALT",
  paste0("var_", samples)
), with = FALSE]

message("Informative SNP/small-indel variants: ", nrow(snp_inf), " / ", nrow(snp))
# -------------------------------------------------------
# load SV
# -------------------------------------------------------
message("Loading SV VCF...")

sv_q <- paste("bcftools query -f '%CHROM\\t%POS\\t%END\\t%ID\\t%INFO/SVTYPE\\t%INFO/SVLEN[\\t%GT]\\n'",shQuote(SV_VCF))

sv <- run_query(sv_q)
setnames(sv, c("chr", "pos", "END", "variant_id", "SVTYPE", "SVLEN", samples))
sv[, pos := as.integer(pos)]
sv[, END := suppressWarnings(as.integer(END))]
sv[, SVLEN := suppressWarnings(as.numeric(SVLEN))]

sv[variant_id == "." | is.na(variant_id), variant_id := paste(chr, pos, SVTYPE, SVLEN, sep = ":")]

sv[, start := pos - 1L]
sv[, end := fifelse(!is.na(END) & END > pos, END, start + pmax(abs(SVLEN), 1, na.rm = TRUE))]
sv[is.na(end) | end <= start, end := start + 1L]

sv[, var_class := SVTYPE]
sv[is.na(var_class) | var_class == ".", var_class := "SV"]

for (s in samples) {
  sv[, paste0("var_", s) := gt_to_present(get(s))]
}

sv[, informative := apply(.SD, 1, is_informative_pattern), .SDcols = paste0("var_", samples)]
sv_inf <- sv[informative == TRUE]

sv_inf <- sv_inf[, c("chr", "start", "end", "variant_id", "var_class", "SVLEN",paste0("var_", samples)), with = FALSE]

# -------------------------------------------------------
# load SVs
# -------------------------------------------------------
message("Informative SVs: ", nrow(sv_inf), " / ", nrow(sv))
message("Informative SNP/small-indel variants: ", nrow(snp_inf), " / ", nrow(snp))

var_dt <- rbindlist(list(snp_inf, sv_inf), use.names = TRUE, fill = TRUE)
var_dt[, variant_group := fifelse(var_class %in% c("SNP", "small_indel"), var_class, "SV")]

fwrite(var_dt, file.path(COMP_DIR, "02_informative_genetic_variants.tsv"), sep = "\t")

find_hits <- function(dmr, var, window_bp = 0L, label = "direct") {
  a <- copy(dmr)
  b <- copy(var)

  a[, query_start := pmax(0L, start - window_bp)]
  a[, query_end := end + window_bp]

  setkey(b, chr, start, end)
  setkey(a, chr, query_start, query_end)

  hits <- foverlaps(a,b,
    by.x = c("chr", "query_start", "query_end"),
    by.y = c("chr", "start", "end"),
    type = "any",
    nomatch = 0L)

  if (nrow(hits) == 0) return(data.table())

  hits[, hit_type := label]
  hits[, distance_bp := fifelse(
    i.start <= end & i.end >= start,
    0L,
    pmin(abs(i.start - end), abs(start - i.end))
  )]

  hits
}

message("Intersecting direct overlaps...")
hits_direct <- find_hits(dmr_dt, var_dt, window_bp = 0L, label = "direct")

message("Intersecting proximal window: ", WINDOW_BP, " bp")
hits_window <- find_hits(dmr_dt, var_dt, window_bp = WINDOW_BP, label = "window")

if (nrow(hits_window) > 0) {
  hits_window[, is_direct_overlap := i.start <= end & i.end >= start]
  hits_prox <- hits_window[is_direct_overlap == FALSE]
} else {
  hits_prox <- data.table()
}

hits_all <- rbindlist(list(hits_direct, hits_prox), use.names = TRUE, fill = TRUE)

score_hits <- function(hits) {
  if (nrow(hits) == 0) return(data.table())

  meth_cols <- paste0("i.meth_", samples)
  var_cols  <- paste0("var_", samples)

  if (!all(meth_cols %in% names(hits))) {
    alt_meth_cols <- paste0("meth_", samples)
    if (all(alt_meth_cols %in% names(hits))) meth_cols <- alt_meth_cols
  }

  out <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
    h <- hits[i]

    meth <- as.integer(unlist(h[, ..meth_cols], use.names = FALSE))
    var  <- as.integer(unlist(h[, ..var_cols], use.names = FALSE))

    sc <- calc_pair_support(meth, var)

    data.table(
      DMRid = h$i.DMRid,
      DMR = h$i.DMR,
      dmr_chr = h$chr,
      dmr_start = h$i.start,
      dmr_end = h$i.end,
      variant_id = h$variant_id,
      variant_class = h$var_class,
      variant_group = h$variant_group,
      var_chr = h$chr,
      var_start = h$start,
      var_end = h$end,
      hit_type = h$hit_type,
      distance_bp = h$distance_bp,
      n_compared = sc$n_compared,
      direct_support = sc$direct_support,
      inverse_support = sc$inverse_support,
      best_support = sc$best_support,
      best_match_frac = sc$best_match_frac,
      support_direction = sc$direction
    )
  }), use.names = TRUE, fill = TRUE)

  out[, passes_support := n_compared >= MIN_SUPPORT &
        best_support >= MIN_SUPPORT &
        best_match_frac >= MIN_MATCH_FRAC]

  out
}

message("Scoring DMR-variant genotype/methylation compatibility...")
pair_scores <- score_hits(hits_all)
#fwrite(pair_scores,file.path(COMP_DIR, sprintf("03_dmr_variant_pair_scores.window_%dbp.tsv", WINDOW_BP)),sep = "\t")
# -------------------------------------------------------
# Test intersect bed
# -------------------------------------------------------
intersect_file <- file.path(COMP_DIR,"06_dmr_genvar.direct_overlap.tsv")
bed_hits <- fread(intersect_file, header = FALSE)

# A = DMR BED: chr start end DMRid DMR
# B = variant BED: chr start end variant_id variant_class REF ALT
setnames(bed_hits, c("dmr_chr", "dmr_start", "dmr_end", "DMRid", "DMR","var_chr", "var_start", "var_end", "variant_id", "variant_class","REF", "ALT"))

var_info <- rbindlist(list(
  snp[, c("variant_id", "var_class", "informative", paste0("var_", samples)), with = FALSE],
  sv[,  c("variant_id", "var_class", "informative", paste0("var_", samples)), with = FALSE]
), use.names = TRUE, fill = TRUE)

setnames(var_info, "var_class", "variant_class")
bed_hits_info <- merge(
  bed_hits,
  var_info,
  by = c("variant_id", "variant_class"),
  all.x = TRUE)

fwrite(bed_hits_info,file.path(COMP_DIR, sprintf("debug_intersectbed_hits_with_informativity.window_%dbp.tsv", WINDOW_BP)),sep = "\t")
