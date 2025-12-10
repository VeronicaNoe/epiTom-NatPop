suppressPackageStartupMessages({
library(data.table)
library(dplyr)
})

setwd("/mnt/disk3/vibanez/metaplot_tomato_group/merged_data")
outDir<-"/mnt/disk3/vibanez/metaplot_tomato_group/"

args <- commandArgs(trailingOnly = TRUE)
#args<-"PIM_CG_ch01"

colNames <- fread(paste0(args[1], "_colNames.tsv"), header = FALSE, sep = '\t', fill = TRUE, na.strings = "NA", nThread = 20)
tmp <- fread(paste0(args[1],".methylation.merged.bed"), 
               sep = '\t', fill = TRUE, na.strings = "NA", nThread = 40)
colnames(tmp) <- c("chr", 'start', 'end', colNames$V1)
avg_data <- tmp %>%
    mutate(avg = rowMeans(select(., all_of(colNames$V1)), na.rm = TRUE)) %>%
    select(chr, start, end, avg)
rm(tmp)
gc()
# Save the averaged data
fwrite(avg_data, paste0(outDir, args[1], "_avg.methylation"), quote = FALSE, row.names = FALSE, col.names = F, sep = "\t")
# Cleanup
