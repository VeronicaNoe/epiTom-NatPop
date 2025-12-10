#########################################################
## run this in Camus
# is a script 
# ~/bin/toRun_check-lambda.sh {}

setwd("/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/sig")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/be_check-GIF/cc_lambda-tables/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/be_check-GIF/cb_plots/"
suppressPackageStartupMessages({
  library(qqman)
  library(data.table)
  library(dplyr)
  library(tidyr)
})
args <- commandArgs(trailingOnly = TRUE)
# args<-"ch12_CG-DMR_9640101"
dmr <- unlist(strsplit(args[1], "_", fixed = TRUE))[2]
input_file <- paste0(args[1], ".ps.gz")

# Read data
toPlot <- data.table::fread(input_file, sep = "\t", data.table =T, 
                            fill = TRUE, header = FALSE, na.strings = "NA", nThread = 20)
colnames(toPlot)<-c('SNPID','beta','seBeta','P')
toPlot <- toPlot %>%
  mutate(SNP=SNPID)%>%
  separate(SNPID, into = c("CHR", "BP"), sep = ":", convert = TRUE)

# Calculate nAssociation
nAss <- sum(toPlot$P <= 1.27e-07)
# Calculate lambda
chisq <- qchisq(1 - toPlot$P, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
inflationFactor <- paste("Inflation Factor:", round(lambda, 5))
# Create data.table
gif <- data.table(dmr = dmr, nAssociation = nAss, lambda = lambda)
# Write to file
out_file <- paste0(outDir, args[1], "_lambda.tsv")
data.table::fwrite(gif, file = out_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
### do the plots
# Preprocess data
toPlot$CHR<-as.numeric(toPlot$CHR)
toPlot$BP<-as.numeric(toPlot$BP)
 # make the plot
png(paste0(outPlot, args[1],"_manhattanPlot.png"))
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
# Manhattan plot
suppressWarnings({
manhattan(toPlot,  suggestiveline = -log10(1.27e-07)) #0.05/395041
})
# QQ plot with inflation factor label
qq(toPlot$P)
# Get the coordinates for the middle of the plotting area
middle_x <- par("pin")[1] - par("pin")[1] # width
middle_y <- (par("pin")[2] + par("pin")[2])*2  # height
text(x = middle_x, y = middle_y, label = inflationFactor, pos = 4,
col = "red", cex = 1.0)
dev.off()
