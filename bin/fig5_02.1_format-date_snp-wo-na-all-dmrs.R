suppressPackageStartupMessages({
  library(dplyr)
})
set.seed(123)  # To make the random split reproducible

metVCFdir<-"/mnt/disk2/vibanez/06_get-meth-vcf/ac_vcf-metabolome/"
snpVCFdir<-"/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/ab_vcf-metabolome/"
leafMETAdir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype/raw/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/ba_formated-data/"

##############################################################################
########################## Load methylation status
##############################################################################
dmr<-data.table::fread(paste0(metVCFdir,'DMR_general_leaf-metabolome_LD.vcf'),
                              skip="#CHROM", sep = '\t', data.table = FALSE,
                              fill = TRUE, na.string=c("NA"), nThread = 20)
tmpCol <- do.call(rbind, strsplit(dmr$ID, ":"))
tmp <- data.frame(tmpCol, stringsAsFactors = FALSE)
dmr$ID<-paste0(gsub('SL2.50','',tmp$X1),'_',tmp$X3,'_',tmp$X2)
rownames(dmr)<-dmr$ID

infoCol<-c('#CHROM', 'POS', 'ID', 'REF','ALT', 'QUAL', 'FILTER','INFO','FORMAT')
col2keep<-setdiff(colnames(dmr), infoCol)
samplesNames<-gsub('_1','', col2keep)
samplesNames<-gsub('S_','S-', samplesNames)
colnames(dmr)<-c(infoCol, samplesNames)

dmr<-dmr[,c('ID',samplesNames)]
head(dmr)
dmr<- dmr %>%
  mutate(across(-ID, ~case_when(
    . == '1/1' ~ 2,  # Replace '1/1' with 2
    . == '0/0' ~ 0,  # Replace '0/0' with 0
    . == './.' ~ NA_real_,  # Replace './.' with NA
    TRUE ~ as.numeric(.) )))%>%
  select(-ID)
dmr[1:5,1:6]
dmr_clean <- dmr[rowSums(is.na(dmr)) == 0,]
dmr_clean[1:5,1:6]
dim(dmr_clean)
#84277
##############################################################################
########################## Load SNP data
##############################################################################
snps<-data.table::fread(paste0(snpVCFdir,'SNP_general_leaf-metabolome_LD.vcf'),
                              skip="#CHROM", sep = '\t', data.table = FALSE,
                              fill = TRUE, na.string=c("NA"), nThread = 20)

tmpCol <- do.call(rbind, strsplit(snps$ID, ":"))
tmp <- data.frame(tmpCol, stringsAsFactors = FALSE)
snps$ID<-paste0(gsub('SL2.50','',tmp$X1),'_',tmp$X3,'_',tmp$X2)
rownames(snps)<-snps$ID

infoCol<-c('#CHROM', 'POS', 'ID', 'REF','ALT', 'QUAL', 'FILTER','INFO','FORMAT')
col2keep<-setdiff(colnames(snps), infoCol)
samplesNames<-gsub('_1','', col2keep)
samplesNames<-gsub('S_','S-', samplesNames)
colnames(snps)<-c(infoCol, samplesNames)
snps<-snps[,c('ID',samplesNames)]
snps[1:6,1:6]

snps <- snps %>%
  mutate(across(-ID, ~case_when(
    . == '1/1' ~ 2, 
    . == '0/1' | . == '1/0'   ~ 1,
    . == '0/0' ~ 0,  
    . == './.' ~ NA_real_,  
    TRUE ~ as.numeric(.) )))%>%
  select(-ID)
snps[1:5,1:6]
snps_clean <- snps[rowSums(is.na(snps)) == 0,] # remove NA
snps_clean[1:5,1:6]
dim(snps_clean)
#1311
gc()
##############################################################################
######################### load metabolite data
##############################################################################
metaboliteData<- data.table::fread(paste0(leafMETAdir,"leaf.metabolites.tsv"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
head(metaboliteData)
metaboliteData$Accessions<-gsub('_1','', metaboliteData$Accessions)
metaboliteData$Accessions<-gsub('S_','S-', metaboliteData$Accessions)
rownames(metaboliteData)<-metaboliteData$Accessions
metaboliteData<-metaboliteData %>%
  select(-Accessions)
metaboliteData[1:7,1:4]

##############################################################################
######################### common samples 
##############################################################################
toKeep<-intersect(rownames(metaboliteData), colnames(dmr))
metaboliteData<-metaboliteData[toKeep,]

# transpose data
dmr_clean<-dmr_clean[,toKeep]
dmr_clean$Marker <- rownames(dmr_clean)
dmr_clean <- data.table::transpose(dmr_clean, keep.names = "Samples", make.names = "Marker")

snps_clean<-snps_clean[,toKeep]
snps_clean$Marker <- rownames(snps_clean)
snps_clean <- data.table::transpose(snps_clean, keep.names = "Samples", make.names = "Marker")
gc()

snps_clean[1:6,1:6]
##############################################################################
######################### save formated data
##############################################################################
# DMR
data.table::fwrite(dmr_clean, paste0(outDir,"DMR_leaf-metabolite.tsv"), quote=F,
                   row.names=T, col.names = T,sep="\t")
# SNPs
data.table::fwrite(snps_clean, paste0(outDir,"SNP_leaf-metabolite.tsv"), quote=F,
                   row.names=T, col.names = T,sep="\t")
# Trais
data.table::fwrite(metaboliteData, paste0(outDir,"Traits_leaf-metabolite.tsv"), quote=F,
                   row.names=T, col.names = T,sep="\t")
