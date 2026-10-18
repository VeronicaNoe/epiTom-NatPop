suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})
set.seed(123)  # To make the random processes reproducible

metVCFdir <- "/mnt/disk2/vibanez/06_get-meth-vcf/ac_vcf-metabolome/"
snpVCFdir <- "/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/ab_vcf-metabolome/"
leafMETAdir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype/raw/"
outDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/ba_formated-data/"

##############################################################################
# Load and Clean Methylation (DMR) Data
##############################################################################
dmr <- data.table::fread(paste0(metVCFdir, 'DMR_general_leaf-metabolome_LD.vcf'),
                         skip = "#CHROM", sep = '\t', data.table = FALSE, fill = TRUE, 
                         na.string = c("NA"), nThread = 20)
tmpCol <- do.call(rbind, strsplit(dmr$ID, ":"))
dmr$ID <- paste0(gsub('SL2.50', '', tmpCol[, 1]), '_', tmpCol[, 3], '_', tmpCol[, 2])
rownames(dmr) <- dmr$ID

infoCol <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
col2keep <- setdiff(colnames(dmr), infoCol)
samplesNames <- gsub('_1', '', col2keep)
samplesNames <- gsub('S_', 'S-', samplesNames)
colnames(dmr) <- c(infoCol, samplesNames)
dmr <- dmr[, c('ID', samplesNames)]
dmr <- dmr %>%
  mutate(across(-ID, ~case_when(
    . == '1/1' ~ 2, 
    . == '0/0' ~ 0,
    . == './.' ~ NA_real_,
    TRUE ~ as.numeric(.)
  ))) %>%
  select(-ID)
dmr_clean <- dmr[rowSums(is.na(dmr)) == 0,]

##############################################################################
# Load and Clean SNP Data
##############################################################################
snps <- data.table::fread(paste0(snpVCFdir, 'SNP_general_leaf-metabolome_LD.vcf'),
                          skip = "#CHROM", sep = '\t', data.table = FALSE, fill = TRUE,
                          na.string = c("NA"), nThread = 20)
tmpCol <- do.call(rbind, strsplit(snps$ID, ":"))

tmp <- data.frame(tmpCol, stringsAsFactors = FALSE)
snps$ID<-paste0(gsub('SL2.50','',tmp$X1),'_',tmp$X3,'_',tmp$X2)
rownames(snps)<-snps$ID

infoCol <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
col2keep <- setdiff(colnames(snps), infoCol)
samplesNames <- gsub('_1', '', col2keep)
samplesNames <- gsub('S_', 'S-', samplesNames)
colnames(snps) <- c(infoCol, samplesNames)

snps <- snps[, c('ID', samplesNames)]
snps[1:6,1:6]

snps <- snps %>%
  mutate(across(-ID, ~case_when(
    . == '1/1' ~ 2, 
    . == '0/1' | . == '1/0' ~ 1,
    . == '0/0' ~ 0,
    . == './.' ~ NA_real_,
    TRUE ~ as.numeric(.)
  ))) %>%
  select(-ID)
snps_clean <- snps[rowSums(is.na(snps)) == 0,] # remove NA
snps_clean[1:5,1:6]
gc()

##############################################################################
# Load Metabolite Data
##############################################################################
metaboliteData <- data.table::fread(paste0(leafMETAdir, "leaf.metabolites.tsv"), 
                                    sep = '\t', data.table = FALSE, fill = TRUE, 
                                    na.string = c("NA"), nThread = 20)
metaboliteData$Accessions <- gsub('_1', '', metaboliteData$Accessions)
metaboliteData$Accessions <- gsub('S_', 'S-', metaboliteData$Accessions)
rownames(metaboliteData) <- metaboliteData$Accessions
metaboliteData <- metaboliteData %>%
  select(-Accessions)

##############################################################################
# Feature Selection
##############################################################################
# # Calculate variance per marker (row) across all samples (columns)
# marker_variances <- apply(dmr_clean, 1, var, na.rm = TRUE)
# # Convert to a data frame for easier viewing
# variance_df <- data.frame(Marker = rownames(dmr_clean), Variance = marker_variances)
# # View the top markers by variance
# variance_df <- variance_df[order(variance_df$Variance, decreasing = TRUE), ]
# # Display the top 10 markers with the highest variance
# head(variance_df, 10)
# tail(variance_df, 10)

# Determine the minimum number of features to balance SNP and DMR datasets
n_features <- min(nrow(snps_clean), nrow(dmr_clean))

# Sort SNPs and DMRs by variance and select the top markers
top_snps <- snps_clean[order(apply(snps_clean, 1, var, na.rm = TRUE), decreasing = TRUE), ][1:n_features, ]
top_dmr <- dmr_clean[order(apply(dmr_clean, 1, var, na.rm = TRUE), decreasing = TRUE), ][1:n_features, ]

##############################################################################
# Save Formatted Data
##############################################################################
# Transpose data for prediction models
toKeep <- intersect(rownames(metaboliteData), colnames(top_snps))
metaboliteData <- metaboliteData[toKeep, ]

top_snps <- top_snps[, toKeep]
top_snps$Marker <- rownames(top_snps)
top_snps <- data.table::transpose(top_snps, keep.names = "Samples", make.names = "Marker")

top_dmr <- top_dmr[, toKeep]
top_dmr$Marker <- rownames(top_dmr)
top_dmr <- data.table::transpose(top_dmr, keep.names = "Samples", make.names = "Marker")

# Save SNP and DMR data
data.table::fwrite(top_snps, paste0(outDir, "SNP_leaf-metabolite.tsv"), quote = FALSE, 
                   row.names = TRUE, col.names = TRUE, sep = "\t")
data.table::fwrite(top_dmr, paste0(outDir, "DMR_leaf-metabolite.tsv"), quote = FALSE, 
                   row.names = TRUE, col.names = TRUE, sep = "\t")
data.table::fwrite(metaboliteData, paste0(outDir, "Traits_leaf-metabolite.tsv"), quote = FALSE, 
                   row.names = TRUE, col.names = TRUE, sep = "\t")
