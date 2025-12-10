suppressPackageStartupMessages({
  library(progress)
  library(data.table)
  library(dplyr)
  library(purrr)
  library(future)
  library(furrr)
})

setwd("/mnt/disk3/vibanez/DMR-GWAS/data")
outDir <- "/mnt/disk3/vibanez/DMR-GWAS/results/"
args <- commandArgs(trailingOnly = TRUE)
#args <- "ch01_C-DMR"

# Read data
input <- list.files(pattern = args[1], full.names = TRUE)
cat("- Step 1: Loading files\n")
plan(multisession, workers = 50)
# Function to process each file
process_file <- function(file_path) {
  dmr <- unlist(strsplit(basename(file_path), ".", fixed = TRUE))[1]
  out <- fread(file_path, sep = '\t', header = FALSE, fill = TRUE, na.strings = "", select = c(1, 4))
  out$dmr <- dmr
  return(out)
}
# Read files in parallel
dt <- future_map_dfr(input, process_file, .progress = TRUE)
gc()
#
cat("- Step 2: Merging DMRs\n")
# Define a function to combine p-values using Fisher's method
combine_pvalues_fisher <- function(pvalues) {
  combined_chisq = -2 * sum(log(pvalues))
  df = 2 * length(pvalues)
  combined_pvalue = pchisq(combined_chisq, df, lower.tail = FALSE)
  min_pvalue = min(pvalues)
  return(list(min_pvalue = min_pvalue,
              combined_pvalue = combined_pvalue, 
              combined_chisq = combined_chisq))
}
# Define a function to merge positions within a specified distance
#setorder(dt,CHR, POS)
setorder(dt,V1,dmr)
# Define the number of chunks
num_chunks <- 100
# Split the data into chunks
chunk_size <- ceiling(nrow(dt) / num_chunks)
chunks <- split(dt, ceiling(seq_len(nrow(dt)) / chunk_size))
# Function to process each chunk
process_chunk <- function(chunk) {
  chunk %>%
    group_by(V1) %>%
    summarize(
	SNPs = first(V1),
	nDMR=n(),
	DMRminPvalue = dmr[which.min(V4)],
	combined_values = list(combine_pvalues_fisher(V4))) %>%
    mutate(
      combined_chisq = map_dbl(combined_values, "combined_chisq"),
      min_p_value = map_dbl(combined_values, "min_pvalue"),
      combined_pvalue = map_dbl(combined_values, "combined_pvalue")
    ) %>%
    select(-combined_values)
}
# Parallel processing of chunks
plan(multisession, workers = 40)
result_chunks <- future_map(chunks, process_chunk)
# Combine the results from all chunks
merged_positions <- bind_rows(result_chunks)

#######
merged_positions<-merged_positions %>%
  mutate(SNPs = V1) %>%
  tidyr::separate(V1, into = c("CHR", "POS"), sep = ":") %>%
  mutate(CHR = as.integer(CHR), POS = as.integer(POS))
# Apply the function to merge positions within 500000bp
data.table::fwrite(merged_positions, file= paste0(outDir,args[1],'.merged.ps'), quote=F,row.names=F,col.names = F, sep="\t")
