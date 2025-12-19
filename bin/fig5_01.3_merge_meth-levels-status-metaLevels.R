# analyze the significant associations
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
  })
  basedir<-"10_data-analysis/Fig5/aa_GWAS-metabolome"
  outDir<-"10_data-analysis/Fig5/ab_data-analysis/results/"
  outPlot<-"10_data-analysis/Fig5/ab_data-analysis/plots/"
  
  ######################### load methylation levles
  inDir<-"05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation/"
  methLevels<-data.table::fread(paste0(inDir,'DMR_methylationLevels'), sep = '\t', 
                                data.table = FALSE,
                                fill = TRUE, na.string=c("NA"), nThread = 20)
  colnames(methLevels)<-gsub('_leaf','',colnames(methLevels))
  colnames(methLevels)<-gsub("LA2838A","S-lycLA2838A",colnames(methLevels))
  colnames(methLevels)<-gsub("TS-677","S-lycLYC3153",colnames(methLevels))
  rownames(methLevels)<-methLevels$DMRs
  head(methLevels)
  ########################## Load methylation status
  methStatus<-data.table::fread("06_get-meth-vcf/ac_vcf-metabolome/DMR_general_leaf-metabolome_LD.vcf", 
                                skip="#CHROM", sep = '\t', data.table = FALSE,
                                fill = TRUE, na.string=c("NA"), nThread = 20)
  tmpCol <- do.call(rbind, strsplit(methStatus$ID, ":"))
  tmp <- data.frame(tmpCol, stringsAsFactors = FALSE)
  methStatus$ID<-paste0(gsub('SL2.50','',tmp$X1),'_',tmp$X3,'_',tmp$X2)
  rownames(methStatus)<-methStatus$ID
  
  infoCol<-c('#CHROM', 'POS', 'ID', 'REF','ALT', 'QUAL', 'FILTER','INFO','FORMAT')
  col2keep<-setdiff(colnames(methStatus), infoCol)
  samplesNames<-gsub('_1','', col2keep)
  samplesNames<-gsub('S_','S-', samplesNames)
  colnames(methStatus)<-c(infoCol, samplesNames)
  methStatus[1:3,1:10]
  
  ######################### load metabolite data
  metaboliteData<- data.table::fread("10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype/raw/leaf.metabolites.tsv", sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  head(metaboliteData)
  metaboliteData$Accessions<-gsub('_1','', metaboliteData$Accessions)
  metaboliteData$Accessions<-gsub('S_','S-', metaboliteData$Accessions)
  rownames(metaboliteData)<-metaboliteData$Accessions
  metaboliteData <- metaboliteData %>%
    mutate(across(-Accessions, log10))
  metaboliteData[1:7,1:4]
  
  ########################### load sig positions
  assoMM<- data.table::fread(paste0(outDir,'01.1_QTLs.tsv'),
                             sep = '\t', data.table = TRUE, fill = TRUE,header = T,
                             na.string=c("NA"), nThread = 20)
  assoMM_DMR <- assoMM %>%
    filter(mm=="DMR" & kinship=="DMR-SNP" & index=="aBN" & resultType=="sig") %>%
    separate(qtl, into = c("ch", "pos", "DMR"), sep = ":") %>%
    mutate(QTL = paste0(ch, "_", DMR, "_", pos)) %>%
    select(-ch, -pos, -DMR)  
  
  #################################################################################                   
  ########################### merge all
  #################################################################################
  mStatus_filtered<-methStatus[unique(assoMM_DMR$QTL),c('ID',samplesNames)]
  mStatus_filtered<-mStatus_filtered%>%
    tidyr::pivot_longer(cols=-c("ID"),
                        names_to = "Accessions",  # New column for metabolite names
                        values_to = "methStatus")
  mStatus_filtered
  mStatus_filtered <- mStatus_filtered %>%
    mutate(
      methStatus = case_when(
        methStatus == '1/1' ~ 2,       
        methStatus == '0/0' ~ 0,       
        methStatus == './.' ~ NA_real_,
        TRUE ~ 3
      )
    )
  #CHECK if there is a 3, there are values different to 1/1,0/0 or ./.
  unique(mStatus_filtered$methStatus)
  
  mStatus_filtered <- mStatus_filtered %>%
    inner_join(assoMM_DMR, by = c("ID" = "QTL"), relationship = "many-to-many")
  
  ###########################################################################################
  mLevels_filtered<-methLevels[unique(assoMM_DMR$QTL),c('DMRs',samplesNames)]
  mLevels_filtered<-mLevels_filtered%>%
    tidyr::pivot_longer(cols=-c("DMRs"),
                        names_to = "Accessions",  # New column for metabolite names
                        values_to = "methLevels")
  
  mLevels_filtered <- mLevels_filtered %>%
    inner_join(assoMM_DMR, by = c("DMRs" = "QTL"),relationship = "many-to-many")
  
  #####################################
  metaLevels<-metaboliteData[samplesNames,c("Accessions",unique(assoMM_DMR$metabolite))]
  metaLevels<-metaLevels%>%
    tidyr::pivot_longer(cols=-c("Accessions"),
                        names_to = "Metabolite",  # New column for metabolite names
                        values_to = "metaLevel")
  metaLevels <- metaLevels %>%
    inner_join(assoMM_DMR, by = c("Metabolite" = "metabolite"),relationship = "many-to-many")
  
  
  mStatus_filtered <- as.data.table(mStatus_filtered)
  mLevels_filtered <- as.data.table(mLevels_filtered)
  metaLevels <- as.data.table(metaLevels)
  
  #### final merge
  merged_data <- mStatus_filtered[
    mLevels_filtered[, .(DMRs, Accessions, metabolite, methLevels)], 
    on = .(ID = DMRs, Accessions, metabolite)
  ]
  
  metaLevels_unique <- metaLevels %>%
    distinct(ID = QTL, Accessions, Metabolite, metaLevel)
  
  merged_data <- merged_data[
    metaLevels_unique,
    on = .(ID, Accessions, metabolite = Metabolite),
    metaLevel := i.metaLevel  # Add the metaLevel column from metaLevels
  ]
  
  data.table::fwrite(merged_data, file= paste0(outDir,"02.0_merged-mLevels-mStatus-metaLevels.tsv"), quote=F,
                     row.names=F,col.names = T,sep="\t")
  QTLcoordinates <-  merged_data %>%
    distinct(ID) %>%                       
    separate(ID, into = c("chr", NA, "start"), sep = "_", remove = FALSE) %>% 
    mutate(
	chr = gsub("^ch", "", chr),
	start = as.numeric(start),  # Ensure `start` is numeric
	end = start + 99            # Calculate `end`
    ) %>%
    select(chr, start, end, ID)   # Reorder columns
  data.table::fwrite(QTLcoordinates, file= paste0(outDir,"02.1_QTL_coordinates.bed"), quote=F,
                     row.names=F,col.names = F,sep="\t")
