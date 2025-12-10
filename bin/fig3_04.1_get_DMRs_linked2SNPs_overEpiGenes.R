#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
})
################################################################################
############ A - linked SNPs per DMR
################################################################################
{
  #results from ~/bin/check_top5chisq.sh
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/associated-SNP-over-epigenes")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results"
input<-list.files(pattern = '.snpsOverEpiGenes', full.names = FALSE)
# Read and process filess
empty_files <- list()
# Process each file
dt <- lapply(input, function(input) {
  dmrPos <- unlist(strsplit(basename(input), ".", fixed = TRUE))[1]
  dmr <- unlist(strsplit(dmrPos, "_", fixed = TRUE))[2]
  # Check if the file is empty before reading
  if (file.size(input) == 0) {
    empty_files <<- c(empty_files, dmr)
    return(NULL)
  }
  out <- fread(input, sep = '\t', header = FALSE, fill = TRUE, na.strings = "")
  out$dmrPos <- dmrPos
  out$dmr <- dmr
  return(out)
})
dt <- Filter(Negate(is.null), dt)
dt <- rbindlist(dt, fill = TRUE)
colnames(dt)<-c('SNP', 'beta','betaSD','pvalue','DMRpos','DMR')
setorder(dt, SNP)

# Check the combined values are similar using Fisher's method
{
# combine_pvalues <- function(pvalues) {
#   chisq_stat <- -2 * sum(log(pvalues))
#   df <- 2 * length(pvalues)  # degrees of freedom
#   p_combined <- pchisq(chisq_stat, df, lower.tail = FALSE)
#   return(p_combined)
# }
# dt_combined <- dt[, .(
#   combined_pvalue = combine_pvalues(pvalue)), by = .(SNP, DMR)]
##
}
onlyLinkedDMRs<- dt %>%
  filter(pvalue <= 1.27e-06)

summary_onlyLinkedDMRs<-dt %>%
  filter(pvalue <= 1.27e-06) %>%
  group_by(SNP)%>%
  dplyr::summarise(avgBeta=mean(beta),
                   count=n())%>%
  arrange(count)
data.table::fwrite(summary_onlyLinkedDMRs,
                   file=paste0(outDir,'04.1_DMRs_linked_summary.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")

data.table::fwrite(as.data.table(unique(onlyLinkedDMRs$DMRpos)),
                   file=paste0(outDir,'04.2_DMRs_linked_list.tsv'), quote=F,
                   row.names=F,col.names = F, sep="\t")

data.table::fwrite(onlyLinkedDMRs,
                   file=paste0(outDir,'04.3_DMRs_linked.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")
}
################################################################################
############ B - meth levels of DMRs
################################################################################
{#try with binary
  input <- fread(paste0(outDir,"04.2_DMRs_linked_list.tsv"), 
                 sep = '\t', header = FALSE, fill = TRUE, na.strings = "")
  dmrDir<-"mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles/"
  #dmrDir<-"/mnt/disk3/vibanez/DMR-GWAS/results/linked_DMRs/"
  input<-list.files(path=dmrDir, pattern = 'general', full.names = FALSE)
  #toKeep<-unique(onlyLinkedDMRs$DMRpos)
  
  linkedDMRs<-c()
  for( i in input){
    dmr<-unlist(strsplit(i, ".", fixed = TRUE))[1]
    dmr<-unlist(strsplit(dmr, "_", fixed = TRUE))[2]
    tmp <- fread(paste0(dmrDir,i), sep = "\t", skip = 1,header = T, fill = TRUE, na.strings = "")
    tmp$DMR<-paste0(gsub('SL2.50','',tmp$`#CHROM`),"_",dmr,"_", tmp$POS)
    tmp<- tmp %>%
      filter(DMR %in% input)
    linkedDMRs<-rbind.data.frame(linkedDMRs, tmp)
  }
  linkedDMRs
  sampleCols<-grep('leaf',colnames(linkedDMRs),value=T)
  infoCol<-setdiff(colnames(linkedDMRs),sampleCols)
  linkedDMRs<-linkedDMRs[, c(infoCol, sampleCols), with = FALSE]
  sampleTab<-data.table::fread("/mnt/disk3/vibanez/DMR-GWAS/results/00_sampleTab.tsv",
                               sep = '\t', data.table = F, fill = TRUE, na.string=c("NA"), nThread = 16)
  sampleTab<-subset(sampleTab, Organ=="leaf")
  rownames(sampleTab)<-sampleTab$MethName
  sampleTab<-sampleTab[sampleCols,]
  snpNames<-sampleTab$SNPnames
  colnames(linkedDMRs)<-c(infoCol, snpNames)

  # from wide to long
  linkedDMRs <- melt(linkedDMRs,
                     id.vars = "DMR",
                     measure.vars = snpNames,
                     variable.name = "Sample",
                     value.name = "methStatus")
  linkedDMRs[, binaryMethStatus := fifelse(methStatus == "0/0", 0,
                                                fifelse(methStatus == "1/1", 2,
                                                        fifelse(methStatus == "./.", NA_real_, NA_real_)))]
  colnames(linkedDMRs)<-c("DMRpos","Sample","MethylationStatus","binary_MethylationStatus")
  data.table::fwrite(linkedDMRs,
                     file=paste0(outDir,'04.4_methStatus_linked_DMRs.tsv'), quote=F,
                     row.names=F,col.names = T, sep="\t")
}

{
  dmrDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bb_phenotype/"
linkedDMRs<-c()
for( i in input){
  dmr<-unlist(strsplit(i, ".", fixed = TRUE))[1]
  tmp <- fread(paste0(dmrDir,i,".phenotype"), sep = "\t", header = F, fill = TRUE, na.strings = "")
  tmp$V1<-dmr
  linkedDMRs<-rbind.data.frame(linkedDMRs, tmp)
}
linkedDMRs
colnames(linkedDMRs)<-c("DMRpos","Sample","MethylationLevels")

data.table::fwrite(linkedDMRs,
                   file=paste0(outDir,'04.5_methLevels_linked_DMRs.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")
}
################################################################################
############ C - genotypes and snp effect over epiGenes
################################################################################
{
snpEffdir<-"/mnt/disk3/vibanez/DMR-GWAS/results/snpEff_genotypes/"
onlyLinkedDMRs[, c("CHR", "pos") := tstrsplit(SNP, ":", fixed = TRUE)]
onlyLinkedDMRs$CHR<-paste0('SL2.50ch',onlyLinkedDMRs$CHR)
onlyLinkedDMRs$pos <- as.integer(onlyLinkedDMRs$pos)

input<-list.files(path=snpEffdir,pattern = '.snpEff_genes_genotypes.txt', full.names = FALSE)
input<-input[sapply(input, function(file) {
  any(sapply(unique(onlyLinkedDMRs$CHR), function(chr) grepl(chr, file)))
})]

snpEffdata<-c()
pattern <- "(?<=\\b[A-Z]/[A-Z])|(?<=\\./\\.)"
for( i in input){
  tmp <- fread(paste0(snpEffdir,i), sep = "\t", header = F, fill = TRUE, na.strings = "")
  filtered_tmp <- merge(tmp, onlyLinkedDMRs, by.x = c("V1", "V2"), by.y = c("CHR", "pos"))
  # get snp effect data
  snpEff_info <- lapply(strsplit(filtered_tmp$V3, "\\|"), function(x) {
    length(x) <- 10  # Adjust this length to the correct number of expected fields
    return(x)
  })
  snpEff_info <- rbindlist(lapply(snpEff_info, function(x) as.list(x)))
  snpEff_info<-snpEff_info[,c(1:4,9,10)]
  setnames(snpEff_info, c("Allele", "Effect", "Impact", "Gene_Name", "nose",
                            "Position_AA_Change"))
  filtered_tmp[, V3 := NULL]

  # now genotypes
  filtered_tmp[, geno_pairs := strsplit(V4, pattern, perl=TRUE)]
  filtered_tmp[, V4 := NULL]
  allInfo<-cbind.data.frame(filtered_tmp, snpEff_info)
  setDT(allInfo)
  genotype_data <- allInfo[, .(SNP, DMRpos, DMR, beta,betaSD, pvalue,
                               Allele, Effect, Impact,Gene_Name, nose,
                               Position_AA_Change, Geno_Pair = unlist(geno_pairs))]
  genotype_data[, c("Sample", "Genotype") := tstrsplit(Geno_Pair, "=", fixed = TRUE)]
  genotype_data[, Geno_Pair := NULL]

  snpEffdata<-rbind.data.frame(snpEffdata, genotype_data)
}
snpEffdata
data.table::fwrite(snpEffdata,file=paste0(outDir,'04.6_SNPs_genotype-effects.tsv'),
                   quote=F, row.names=F,col.names = T, sep="\t")
}
################################################################################
############ D - plot MethLevels & Genotype per Epigene
################################################################################
rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
})
outDir<-"/mnt/disk3/vibanez/DMR-GWAS/results/"
linkedDMRs_levels<- fread(paste0(outDir,"04.5_methLevels_linked_DMRs.tsv"), sep = "\t", header = T,
                          fill = TRUE, na.strings = "")
linkedDMRs_status<- fread(paste0(outDir,"04.4_methStatus_linked_DMRs.tsv"),
                          sep = "\t", header = T, fill = TRUE, na.strings = "")

snpEffdata<- fread(paste0(outDir,"04.6_SNPs_genotype-effects.tsv"), sep = "\t", header = T, fill = TRUE, na.strings = "")

## methLevels
{
toKeep<-intersect(unique(linkedDMRs_levels$Sample), unique(snpEffdata$Sample))
 snpEffdata<-snpEffdata %>%
   filter(Sample %in% toKeep)
 linkedDMRs_levels<-linkedDMRs_levels %>%
   filter(Sample %in% toKeep)

 methGenotype <- merge(snpEffdata, linkedDMRs_levels, by = c("DMRpos","Sample"), all.x = TRUE)
 methGenotype$MethylationLevels <- as.numeric(methGenotype$MethylationLevels)
 methGenotype$Gene_SNP_DMR <- paste(methGenotype$Gene_Name, methGenotype$SNP, methGenotype$DMRpos, sep = "_")
 methGenotype$Gene_SNP <- paste(methGenotype$Gene_Name, methGenotype$SNP, sep = "_")

geneName<- fread(paste0(outDir,"04.7_SNPs_over_epiGenes_annotated.tsv"),
                           sep = "\t", header = T, fill = TRUE, na.strings = "")
merged_data <- merge(geneName, methGenotype, by.x = "geneName", by.y = "Gene_Name", all = TRUE, allow.cartesian = TRUE)

merged_data<-merged_data %>%
   filter(is.na(Name)!=1)
# Generate the violinplot
for(i in unique(merged_data$Gene_SNP)){
   tmp<-merged_data %>%
     filter(Gene_SNP==i & Genotype!="./.")
   snp<-unique(tmp$SNP)
   dmr<-unique(tmp$DMRpos)
   gene<-unique(tmp$Gene_Name)
   ggplot(tmp, aes(x = Genotype, y = MethylationLevels)) +
     geom_violin(trim = FALSE, fill = "lightblue", color = "black") +  # Violin plot
     geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "darkblue") +  # Overlay jittered points
     labs(
       x = paste0(snp, " genotype over ", gene ),
       y = paste0(dmr, " Methylation Levels"),
       title = paste0("DMR: ", dmr," - SNP: ", snp, " - Gene: ", gene)
     ) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
   ggsave(paste0(outPlot,"/tmp/04.2_methLevels_violinPlot_",i,".pdf"), width = 12, height = 8)
 }
 #for(i in unique(methGenotype$Gene_SNP_DMR)){
#   tmp<-methGenotype %>%
#     filter(Gene_SNP_DMR==i & Genotype!="./.")
#   snp<-unique(tmp$SNP)
#   dmr<-unique(tmp$DMRpos)
#   gene<-unique(tmp$Gene_Name)
#   ggplot(tmp, aes(x = Genotype, y = MethylationLevels)) +
#     geom_violin(trim = FALSE, fill = "lightblue", color = "black") +  # Violin plot
#     geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "darkblue") +  # Overlay jittered points
#     labs(
#       x = paste0(snp, " genotype over ", gene ),
#       y = paste0(dmr, " Methylation Levels"),
#       title = paste0("DMR: ", dmr," - SNP: ", snp, " - Gene: ", gene)
#     ) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
#   ggsave(paste0("plots/03.2_methLevels_violinPlot_",i,".pdf"), width = 12, height = 8)
# }
}
# status
toKeep<-intersect(unique(linkedDMRs_status$Sample), unique(snpEffdata$Sample))
snpEffdata<-snpEffdata %>%
  filter(Sample %in% toKeep)
linkedDMRs_status<-linkedDMRs_status %>%
  filter(Sample %in% toKeep)

methGenotype <- merge(snpEffdata, linkedDMRs_status, by = c("DMRpos","Sample"), all.x = TRUE)
methGenotype$binary_MethylationStatus <- as.numeric(methGenotype$binary_MethylationStatus)
methGenotype$Gene_SNP_DMR <- paste(methGenotype$Gene_Name, methGenotype$SNP, methGenotype$DMRpos, sep = "_")
methGenotype$Gene_SNP <- paste(methGenotype$Gene_Name, methGenotype$SNP, sep = "_")
# Generate the violinplot
for(i in unique(methGenotype$Gene_SNP)){
  tmp<-methGenotype %>%
    filter(Gene_SNP==i & Genotype!="./.")
  snp<-unique(tmp$SNP)
  dmr<-unique(tmp$DMRpos)
  gene<-unique(tmp$Gene_Name)
  ggplot(tmp, aes(x = Genotype, y = binary_MethylationStatus)) +
    geom_violin(trim = FALSE, fill = "lightblue", color = "black") +  # Violin plot
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "darkblue") +  # Overlay jittered points
    labs(
      x = paste0(snp, " genotype over ", gene ),
      y = paste0(dmr, " Methylation Levels"),
      title = paste0("SNP: ", snp, " - Gene: ", gene)
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
  ggsave(paste0(outPlot,"/tmp/04.2_methStatus_violinPlot_",i,".pdf"), width = 12, height = 8)
}
for(i in unique(methGenotype$Gene_SNP_DMR)){
  tmp<-methGenotype %>%
    filter(Gene_SNP_DMR==i & Genotype!="./.")
  snp<-unique(tmp$SNP)
  dmr<-unique(tmp$DMRpos)
  gene<-unique(tmp$Gene_Name)
  ggplot(tmp, aes(x = Genotype, y = binary_MethylationStatus)) +
    geom_violin(trim = FALSE, fill = "lightblue", color = "black") +  # Violin plot
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "darkblue") +  # Overlay jittered points
    labs(
      x = paste0(snp, " genotype over ", gene ),
      y = paste0(dmr, " Methylation Levels"),
      title = paste0("DMR: ", dmr," - SNP: ", snp, " - Gene: ", gene)
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
  ggsave(paste0(outPlot,"/tmp/04.3_methStatus_violinPlot_",i,".pdf"), width = 12, height = 8)
}
}
