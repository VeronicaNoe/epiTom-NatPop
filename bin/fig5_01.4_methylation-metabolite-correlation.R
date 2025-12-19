suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
  })
  
  basedir<-"10_data-analysis/Fig5/aa_GWAS-metabolome"
  outDir<-"10_data-analysis/Fig5/ab_data-analysis/results/"
  outPlot<-"10_data-analysis/Fig5/ab_data-analysis/plots/"
  df<-data.table::fread(paste0(outDir,"02.0_merged-mLevels-mStatus-metaLevels.tsv"),
                        sep = '\t', data.table = TRUE, 
                        fill = TRUE, na.string=c("NA"), nThread = 20)
  
  meta <- unique(df$metabolite)
  corOut <- list()  # Use a list to store results efficiently
  
  for (m in meta) {
    subTmp <- df %>%
      filter(metabolite == m)
    
    dmrPositions <- unique(subTmp$ID)  # Correct variable usage
    
    for (d in dmrPositions) {
      tmpDF <- subTmp %>%
        filter(ID == d & methStatus != './.')
      
      # Initialize default values
      pearsonMetaMeth <- NA
      pVal <- NA
      pValueWilcoxTest.MetaMeth <- NA
      
      # Correlation test
      if (n_distinct(tmpDF$methStatus) >= 2) {
        corMetaMeth <- cor.test(tmpDF$metaLevel, tmpDF$methLevels)
        pVal <- corMetaMeth$p.value
        pearsonMetaMeth <- as.numeric(corMetaMeth$estimate)
      }
      
      # Wilcoxon test
      if (n_distinct(tmpDF$methStatus) >= 2) {
        pValueWilcoxTest.MetaMeth <- wilcox.test(metaLevel ~ methStatus, 
                                                 data = tmpDF, exact = FALSE)$p.value
      }
      
      # Append results
      corOut[[length(corOut) + 1]] <- data.frame(
        metabolite = m, 
        DMR = d, 
        pearsonMetaMeth = pearsonMetaMeth, 
        pVal = pVal, 
        pValueWilcoxTest.MetaMeth = pValueWilcoxTest.MetaMeth
      )
    }
  }
  
  # Combine results into a single data frame
  corOut <- do.call(rbind, corOut)
  data.table::fwrite(as.data.frame(corOut), 
                   file= paste0(outDir,'03.0_methylation-metabolite_correlation.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")
