#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
suppressPackageStartupMessages({
  library(methylKit)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: script.R <controlAcc> <sampleAcc> <context: CG|CHG|CHH> <contig>")
}

# =========================================================
# Paths
# =========================================================

path2acc <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ag_filter/"
path2results <- "/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ah_methylome-comparison/"

dir.create(path2results, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Arguments
# =========================================================

controlAcc <- args[1]
sampleAcc  <- args[2]
contextArg <- args[3]
contigArg  <- args[4]
#test
# controlAcc <- "G0-C3"
# sampleAcc  <- "G0-C6"
# contextArg <- "CG"
# contigArg  <- "ptg000152l"

CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "20"))
CHUNK <- 1e6

refGen <- "cervil_polished"

cat("----- Control:", controlAcc, "\n")
cat("----- Sample:", sampleAcc, "\n")
cat("----- Context:", contextArg, "\n")
cat("----- Contig:", contigArg, "\n")
cat("----- Cores:", CORES, "\n")

# =========================================================
# Context conversion
# =========================================================
methylContext <- switch(
  contextArg,
  CG  = "CpG",
  CHG = "CHG",
  CHH = "CHH",
  stop("Unknown context. Use CG, CHG, or CHH.")
)

# =========================================================
# Find input files
# Expected filenames should contain:
# sample, context, contig, and end with .filtered.bed
# Example:
# G0-C3_CHH_ptg000001l.filtered.bed
# =========================================================
find_file <- function(acc, contextArg, contigArg) {
  x <- list.files(
    path = path2acc,
    pattern = "\\.filtered\\.bed$",
    full.names = FALSE
  )
  
  x <- x[
    grepl(acc, x, fixed = TRUE) &
      grepl(contextArg, x, fixed = TRUE) &
      grepl(contigArg, x, fixed = TRUE)
  ]
  
  if (length(x) != 1) {
    cat("Files found for:", acc, contextArg, contigArg, "\n")
    print(x)
    stop("Expected exactly one input file.")
  }
  
  file.path(path2acc, x)
}

controlFile <- find_file(controlAcc, contextArg, contigArg)
sampleFile  <- find_file(sampleAcc,  contextArg, contigArg)

inFile <- c(controlFile, sampleFile)

sampleNames <- c(
  paste(controlAcc, contextArg, contigArg, sep = "_"),
  paste(sampleAcc,  contextArg, contigArg, sep = "_")
)

# group 0 = control, group 1 = sample
contr <- c(0, 1)

cat("----- Input files:\n")
print(inFile)

# =========================================================
# Output files
# =========================================================
out_filtered <- paste0(
  path2results,
  controlAcc, "_", sampleAcc, "_",
  contextArg, "_", contigArg,
  "_DMR.methylation.bed"
)

out_all <- paste0(
  path2results,
  controlAcc, "_", sampleAcc, "_",
  contextArg, "_", contigArg,
  "_DMR.wo-filtering"
)

# =========================================================
# Read methylation data
# =========================================================
cat("----- Loading data to methylKit\n")

cov <- methRead(
  as.list(inFile),
  context = methylContext,
  sample.id = as.list(sampleNames),
  treatment = contr,
  assembly = refGen,
  pipeline = "bismarkCytosineReport",
  mincov = 5,
  header = FALSE
)

# =========================================================
# Coverage filtering
# =========================================================

cat("----- Filtering coverage\n")

filtered <- filterByCoverage(
  cov,
  lo.count = 5,
  lo.perc = NULL,
  hi.count = NULL,
  hi.perc = 99.9,
  chunk.size = CHUNK
)

# =========================================================
# Tile first, then unite
# =========================================================

cat("----- Creating 100 bp windows per sample\n")

wind <- tileMethylCounts(
  filtered,
  win.size = 100,
  step.size = 100,
  cov.bases = 5,
  mc.cores = CORES
)

cat("----- Merging window-level data\n")

norm <- unite(
  wind,
  destrand = FALSE,
  min.per.group = 1L
)

n_norm <- nrow(getData(norm))
cat("----- Windows after unite:", n_norm, "\n")

if (n_norm == 0) {
  stop("No windows remained after unite(). Check coverage, context, contig, and input files.")
}

# =========================================================
# Differential methylation
# =========================================================

cat("----- Calculating DMRs\n")

DMR_obj <- calculateDiffMeth(
  norm,
  mc.cores = CORES
)

# =========================================================
# Extract data safely
# =========================================================

cat("----- Extracting methylKit data\n")

dfMeth <- as.data.table(getData(norm))
dfStat <- as.data.table(getData(DMR_obj))

cat("Rows methylation table:", nrow(dfMeth), "\n")
cat("Rows DMR statistics:", nrow(dfStat), "\n")

# =========================================================
# Merge methylation counts and DMR statistics
# =========================================================

merge_cols <- intersect(
  c("chr", "start", "end", "strand"),
  intersect(names(dfMeth), names(dfStat))
)

cat("Merge columns:", paste(merge_cols, collapse = ", "), "\n")

stat_cols <- c(merge_cols, "pvalue", "qvalue", "meth.diff")
stat_cols <- intersect(stat_cols, names(dfStat))

DMR_all <- merge(
  dfMeth,
  dfStat[, ..stat_cols],
  by = merge_cols,
  all = FALSE
)

cat("Rows after merge:", nrow(DMR_all), "\n")

# =========================================================
# Filter significant DMRs
# =========================================================

threshold <- ifelse(contextArg == "CHH", 10, 50)

DMR_filtered <- DMR_all[
  !is.na(qvalue) &
    qvalue <= 0.05 &
    !is.na(meth.diff) &
    (meth.diff <= -threshold | meth.diff >= threshold)
]

cat("Rows after DMR filtering:", nrow(DMR_filtered), "\n")

# =========================================================
# Write outputs
# =========================================================

cat("----- Writing filtered DMRs\n")

fwrite(
  DMR_filtered,
  file = out_filtered,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t",
  nThread = CORES,
  buffMB = 128
)

cat("----- Writing all tested windows\n")

fwrite(
  DMR_all,
  file = out_all,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t",
  nThread = CORES,
  buffMB = 128
)

cat("----- Done\n")
cat("Filtered:", out_filtered, "\n")
cat("All:", out_all, "\n")
