suppressPackageStartupMessages({
	library(data.table)
	library(furrr)
	library(future)
	library(dplyr)
})
setwd("/mnt/disk3/vibanez/DMR-GWAS/data")
outDir <- "/mnt/disk3/vibanez/DMR-GWAS/results/"
args <- commandArgs(trailingOnly = TRUE)
#args <- "C-DMR"
# Define the number of CPUs to use
# List all files
input_files <- list.files(pattern = args[1], full.names = TRUE)
#input_files <- head(input_files)
# Define the number of CPUs to use
num_cpus <- 10  # Adjust this based on your system's capabilities
# Set up parallel processing with the specified number of CPUs
plan(multisession, workers = num_cpus)
# Function to read and annotate a chunk of rows from each file
read_chunk <- function(file_path, chunk_size, chunk_num) {
  dt <- fread(file_path, sep = '\t', header = FALSE, fill = TRUE, na.strings = "",
              skip = (chunk_num - 1) * chunk_size, nrows = chunk_size)
  if (nrow(dt) == 0) return(NULL)
  dt$dmr <- gsub('.ps.gz','',basename(file_path))
  return(dt)
}
# Define the chunk size
chunk_size <- 100000  # Adjust this based on your memory capacity
# Estimate the number of chunks
# Let's assume each file has approximately 1.2 million rows
approx_rows_per_file <- 1300000
max_chunks <- ceiling(approx_rows_per_file / chunk_size)
# Initialize an empty list to store results
results_list <- list()
# Function to combine p-values using Fisher's method
combine_pvalues_fisher <- function(pvalues) {
  combined_chisq = -2 * sum(log(pvalues))
  df = 2 * length(pvalues)
  combined_pvalue = pchisq(combined_chisq, df, lower.tail = FALSE)
  min_pvalue = min(pvalues)
  return(list(min_pvalue = min_pvalue,
              combined_pvalue = combined_pvalue, 
              combined_chisq = combined_chisq))
}
# Process each chunk in parallel
for (chunk_num in 1:max_chunks) {
  cat(sprintf("\n Processing chunk %d/%d\n", chunk_num, max_chunks))
  # Read chunks from all files in parallel
  dt_list <- future_map(input_files, ~read_chunk(.x, chunk_size, chunk_num), .progress = TRUE)
  # Filter out NULL results (empty chunks)
  dt_list <- dt_list[!sapply(dt_list, is.null)]
  # Combine all data.tables by row (V1)
  combined_dt <- rbindlist(dt_list)
  # Split the combined data.table by V1
  split_dt <- split(combined_dt, by = "V1")
  # Apply the function to each split data.table to compute combined p-values
  chunk_results <- lapply(split_dt, function(dt) {
    pvalues <- dt$V4
    combined_values <- combine_pvalues_fisher(pvalues)
    data.table(
      V1 = dt$V1[1],
      combined_chisq = combined_values$combined_chisq,
      combined_pvalue = combined_values$combined_pvalue,
      min_p_value = combined_values$min_pvalue,
      DMRminPvalue = dt$dmr[which.min(dt$V4)],
      nDMR=nrow(dt)
    )
  })
  # Store the chunk results
  results_list[[chunk_num]] <- rbindlist(chunk_results)
}
# Combine all chunk results into a final data.table
final_results <- rbindlist(results_list)
# Write the final results to a file
fwrite(final_results, file = "combined_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
