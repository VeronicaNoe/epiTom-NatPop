suppressPackageStartupMessages({
  library(qqman)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)

indir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/sig"
outPlot <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/be_check-GIF/cb_plots"
outDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/be_check-GIF/cc_lambda-tables"
file <- args[1]
#input <- list.files(path = indir, pattern = ".ps.gz", full.names = TRUE)
#input<-grep("Tm078", input, value=T)
process_file <- function(filepath) {
  file_parts <- strsplit(basename(filepath), '[._]', fixed = FALSE)[[1]]
  metabolite <- file_parts[1]
  tissue <- file_parts[2]
  molmarker <- file_parts[3]
  kin <- file_parts[4]
  index <- file_parts[5]
  outName <- paste0(metabolite, '_', tissue,'_',molmarker, '_', kin, '_', index)
  
  toPlot <- fread(paste0(indir,"/",file,".ps.gz"), sep = '\t', na.string = c("NA"), nThread = 20,
                  header = FALSE, fill = TRUE)
  
  toPlot[, CHR := as.numeric(sub("ch", "", sapply(strsplit(V1, ":"), `[[`, 1)))]
  toPlot[, BP := as.numeric(sapply(strsplit(V1, ":"), `[[`, 2))]
  setnames(toPlot, c("V1", "V4", "CHR", "BP"), c("SNP", "P", "CHR", "BP"))
  
  nAss <- sum(toPlot$P <= 5e-8)
  chisq <- qchisq(1 - toPlot$P, 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  gif <- data.frame(sample= file, metabolite = metabolite, molmarker = molmarker, kinship=kin, index=index,
                    number_QLTs = nAss, lambda = lambda)
  # Save GIF and lambda file
  fwrite(gif, file = file.path(outDir, paste0(outName, "_lambda.tsv")), quote = FALSE,
         row.names = FALSE, col.names = TRUE, sep = "\t")
  
	toManhatann <- toPlot[-log10(P) > 2, .(SNP, CHR, BP, P)] 
  # Generate plots
png(file = file.path(outPlot, paste0(outName, "_manhattan-plot.png")))
  manhattan(toManhatann, main = outName, suggestiveline = FALSE,ylim = c(2, max(-log10(toPlot$P)+2)))
invisible(dev.off())

png(file = file.path(outPlot, paste0(outName, "_qq-plot.png")))
  qq(toPlot$P, main = paste("QQ Plot -", outName))
  mtext(paste("Genomic Inflation Factor (λ) =", round(lambda, 5)),
      side = 1, line = 2, col = "red", cex = 1.0)
invisible(dev.off())

}

# Apply processing to all files
lapply(file, process_file)
