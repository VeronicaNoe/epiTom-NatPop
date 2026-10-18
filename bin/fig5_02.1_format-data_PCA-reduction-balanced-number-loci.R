suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(FactoMineR)  # For PCA
  library(factoextra)  # For PCA visualization and variance extraction
  library(future.apply)
})

set.seed(123)  # To make the random processes reproducible

# Define directories
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

naThreshold<-0.2*ncol(dmr)
dmr_clean <- dmr[rowSums(is.na(dmr)) <= naThreshold,] # remove NA
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
naThreshold<-0.2*ncol(snps)
snps_clean <- snps[rowSums(is.na(snps)) <= naThreshold,] # remove NA
##############################################################################
# Load and Align Metabolite Data
##############################################################################
metaboliteData <- data.table::fread(paste0(leafMETAdir, "leaf.metabolites.tsv"), 
                                    sep = '\t', data.table = FALSE, fill = TRUE, 
                                    na.string = c("NA"), nThread = 20)
metaboliteData$Accessions <- gsub('_1', '', metaboliteData$Accessions)
metaboliteData$Accessions <- gsub('S_', 'S-', metaboliteData$Accessions)
rownames(metaboliteData) <- metaboliteData$Accessions

# Align with PCA-transformed data
toKeep <- intersect(rownames(metaboliteData), colnames(dmr))
metaboliteData <- metaboliteData[toKeep, ]

##############################################################################
# Balance Feature Count with PCA
##############################################################################
# Handle missing data by median imputation
# Subset and transpose
snps_clean <- snps_clean[, toKeep]
snps_clean$Marker <- rownames(snps_clean)
# Transpose for marker-wise operations
snps_clean <- data.table::transpose(snps_clean, keep.names = "Samples", make.names = "Marker")
snps_clean <- as.data.table(snps_clean)
# Perform imputation over markers (now columns) in parallel

plan(multisession, workers = 80)  # Use 80 parallel workers

snps_clean <- as.data.table(future_lapply(snps_clean, function(col) {
  if (any(is.na(col))) {
    median(col, na.rm = TRUE)
  } else {
    col
  }
}))

dmr_clean <- dmr_clean[, toKeep]
dmr_clean$Marker <- rownames(dmr_clean)
dmr_clean <- data.table::transpose(dmr_clean, keep.names = "Samples", make.names = "Marker")
dmr_clean  <- as.data.table(dmr_clean )
dmr_clean  <- as.data.table(future_lapply(dmr_clean , function(col) {
  if (any(is.na(col))) {
    median(col, na.rm = TRUE)
  } else {
    col
  }
}))

# Perform PCA using prcomp for efficiency
pca_dmr <- prcomp(dmr_clean[, !"Samples", with = FALSE], scale. = TRUE, rank. = 150)
pca_snps <- prcomp(snps_clean[, !"Samples", with = FALSE], scale. = TRUE, rank. = 150)

# Calculate variance explained and find components explaining 90%
explained_variance_dmr <- cumsum(pca_dmr$sdev^2 / sum(pca_dmr$sdev^2) * 100)
explained_variance_snps <- cumsum(pca_snps$sdev^2 / sum(pca_snps$sdev^2) * 100)
top_dmr_components <- which(explained_variance_dmr >= 90)[1]
top_snps_components <- which(explained_variance_snps >= 90)[1]

# Filter loci based on contributions > 0
dmr_contributions <- as.data.frame(pca_dmr$rotation[, 1:top_dmr_components])
snps_contributions <- as.data.frame(pca_snps$rotation[, 1:top_snps_components])
dmr_filtered <- dmr_contributions[rowSums(dmr_contributions > 0) > 0, , drop = FALSE]
snps_filtered <- snps_contributions[rowSums(snps_contributions > 0) > 0, , drop = FALSE]

# Balance loci by retaining the minimum number
min_loci <- min(nrow(dmr_filtered), nrow(snps_filtered))
filtered_dmr <- dmr_filtered[1:min_loci, , drop = FALSE]
filtered_snps <- snps_filtered[1:min_loci, , drop = FALSE]

# Subset original datasets
reduced_dmr <- dmr_clean[, c("Samples", rownames(filtered_dmr)), with = FALSE]
reduced_snps <- snps_clean[, c("Samples", rownames(filtered_snps)), with = FALSE]

# Check dimensions
dim(reduced_dmr)
dim(reduced_snps)

# pca_dmr <- PCA(dmr_clean[, !"Samples", with = FALSE], ncp = 150, graph = FALSE)
# pca_snps <- PCA(snps_clean[, !"Samples", with = FALSE], ncp= 150, graph = FALSE)
# 
# # Retain components explaining 90% variance
# explained_variance_dmr <- cumsum(pca_dmr$eig[, 2])
# explained_variance_snps <- cumsum(pca_snps$eig[, 2])
# top_dmr_components <- which(explained_variance_dmr >= 95)[1]
# top_snps_components <- which(explained_variance_snps >= 95)[1]
# 
# # Convert to data.table for efficiency
# dmr_contributions <- as.data.frame(pca_dmr$var$coord[, 1:top_dmr_components]) %>%
#   mutate(DMR = rownames(pca_dmr$var$coord))%>%
#   filter(if_any(starts_with("Dim."), ~ . > 0))
# 
# snps_contributions <- as.data.frame(pca_snps$var$coord[, 1:top_snps_components]) %>%
#   mutate(SNP = rownames(pca_snps$var$coord))%>%
#   filter(if_any(starts_with("Dim."), ~ . > 0))
# 
# # Subset the original DMR dataset
# reduced_dmr <- dmr_clean[, c("Samples", rownames(dmr_contributions)), with = FALSE]
# reduced_snps <- snps_clean[, c("Samples", rownames(snps_contributions)), with = FALSE]
# dim(reduced_snps)
# dim(reduced_dmr)
# # Calculate the minimum number of loci to retain
# min_loci <- min(ncol(reduced_dmr), ncol(reduced_snps))
# 
# #
# # Select the top `min_loci` DMRs and SNPs
# filtered_dmr <- dmr_contributions %>%
#   arrange(desc(if_any(starts_with("Dim."), ~ .))) %>%
#   slice(1:min_loci)
# 
# filtered_snps <- snps_contributions %>%
#   arrange(desc(if_any(starts_with("Dim."), ~ .))) %>%
#   slice(1:min_loci)
# 
# reduced_dmr <- dmr_clean[, c("Samples", rownames(filtered_dmr )), with = FALSE]
# reduced_snps <- snps_clean[, c("Samples", rownames(filtered_snps)), with = FALSE]
# dim(reduced_snps)
# dim(reduced_dmr)
# 119 86223
# 119 86222



##############################################################################
# Save Outputs
##############################################################################
# Save SNP and DMR data
data.table::fwrite(reduced_snps, paste0(outDir, "SNP_leaf-metabolite.tsv"), quote = FALSE, 
                   row.names = TRUE, col.names = TRUE, sep = "\t")
data.table::fwrite(reduced_dmr, paste0(outDir, "DMR_leaf-metabolite.tsv"), quote = FALSE, 
                   row.names = TRUE, col.names = TRUE, sep = "\t")
data.table::fwrite(metaboliteData, paste0(outDir, "Traits_leaf-metabolite.tsv"), quote = FALSE, 
                   row.names = TRUE, col.names = TRUE, sep = "\t")

