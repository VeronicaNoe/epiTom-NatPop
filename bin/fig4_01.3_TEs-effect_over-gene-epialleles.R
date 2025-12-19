suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(broom)
  
  library(ggpubr)
  library(ggrepel)
  
  library(purrr)
  
  
  library("RColorBrewer")
  library("gplots")
  library(viridis)
  
})
setwd("10_data-analysis/Fig4/aa_get-epialleles-over-genes/bd_merged-gene-epialleles")
outDir<-"10_data-analysis/Fig4/results/"
outPlot<-"10_data-analysis/Fig4/plots/"


allData<-data.table::fread(paste0(outDir,'01.06_epiallele_sample.tsv'), sep = '\t',
                           data.table = T, fill = TRUE, check.names=FALSE,
                           na.string=c("NA"), nThread = 10)
########################################################################################
##  D.1 | TE effect over genes: gene-TE case, DMRs over TEs
######################################################################################## 
## Calculate the effect of DMR over TE
{
  dmrWtes<-unique(subset(allData, teMnoverTEs!=0)$geneName)
  dmrOverTe <- allData
  dmrOverTe[, DMRoverTEs := ifelse(geneName %in% dmrWtes, "overTEs", "woTEs")]
  ## Fisher test
  {
    # fisher test 
    # OR= d/c / b/a in theory
    # OR= a*d / b*c in R
    #       FALSE TRUE
    # FALSE   a   c
    # TRUE    b   d
    # Odd= # enfermos q fueron al zoo / # enfermos que NO fueron
    #       # sanos que fueron al zoo / # sanos que No fueron
    exon <- dmrOverTe %>%
      filter(annotation=="exon" & geneType=="gene-TE")
    contingency_table <- table(exon$DMRoverTEs == "overTEs", exon$sampleEpiallele == "teM")
    fisher.test(contingency_table)

    intron <- dmrOverTe %>%
      filter(annotation=="intron" & geneType=="gene-TE")
    contingency_table <- table(intron$DMRoverTEs == "overTEs", intron$sampleEpiallele == "teM")
    fisher.test(contingency_table)

    total <- dmrOverTe %>%
      filter(geneType=="gene-TE")
    contingency_table <- table(total$DMRoverTEs == "overTEs", total$sampleEpiallele == "teM")
    fisher.test(contingency_table)

  }
  ## get freq
  dmrOverTeFreq <- dmrOverTe %>%
    group_by(geneName,geneEpiallele,geneType, annotation, sampleEpiallele,DMRoverTEs) %>%
    dplyr::summarise(frequency = n()/186,
                     count = n())
  ## plots
  ggplot(dmrOverTeFreq, aes(x = frequency, fill = sampleEpiallele, color=sampleEpiallele, alpha=0.7)) +
    geom_histogram(position = "identity", breaks = seq(0, max(dmrOverTeFreq$frequency), by = 0.01)) +
    geom_freqpoly(position = "identity",bins = 100, binwidth=0.01) +
    scale_fill_manual(values = c('#ffca7b', '#820a86', "grey")) +
    scale_color_manual(values = c('#ffca7b', '#820a86', "grey")) +
    scale_x_continuous(breaks = seq(0, max(dmrOverTeFreq$frequency), by = 0.1)) +  
    labs(title = "Histogram of Allele Frequencies",
         x = "Allele Frequency",
         y = "Number of gene per epiallele") +
    facet_grid(.~DMRoverTEs, scale="free") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  ggsave(paste0(outPlot,"03.01.gene-TE_TEeffect_Epiallele_frequency_histogram.pdf"),  width = 20, height = 20, units = "cm")  
  
  ggplot(dmrOverTeFreq, aes(x = DMRoverTEs, y = count, color=sampleEpiallele, fill=sampleEpiallele)) +
    geom_boxplot(width=0.5, na.rm = T,alpha=0.6,lwd = 0.3,outlier.shape = NA)+
    stat_summary(fun = mean, geom = "point", shape = 16, size = 2,
                 color = "darkred", position = position_dodge(width = 0.5)) +
    #stat_compare_means() +   # Add global p-value
    scale_fill_manual(values = c('#ffca7b','#820a86', "grey")) +
    scale_color_manual(values = c( '#ffca7b','#820a86', "grey")) +
    labs(title="Effect of TEs over gene-TEs", x="", y ="Number of teM epiallele per gene")+
    theme_bw() +
    #facet_grid(annotation~., scale="free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    theme(panel.grid.major = element_line(color = "gray90", linetype = "dashed"))
  ggsave(paste0(outPlot,"03.02.gene-TE_TEeffect_Epiallele_boxplot.pdf"),  width = 20, height = 20, units = "cm")
  
  dmrOverTeBins <- dmrOverTeFreq %>%
    mutate(frequency_bin = cut(frequency, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE))
  # Count the occurrences within each bin
  frequency_counts <- dmrOverTeBins %>%
    #filter(sampleEpiallele!="gbM") %>%
    group_by(geneType, geneEpiallele, sampleEpiallele,annotation, DMRoverTEs, frequency_bin) %>%
    dplyr::summarize(count = n())
  
  ## Get the boxplots of frequencies
  ggplot(frequency_counts, aes(x = frequency_bin, y = count, color=sampleEpiallele, 
                               fill=sampleEpiallele)) +
    geom_boxplot(width=0.5, na.rm = T,alpha=0.6,lwd = 0.3,outlier.shape = NA)+
    stat_summary(fun = mean, geom = "point", shape = 16, size = 2, 
                 color = "darkred", position = position_dodge(width = 0.5)) + 
    #stat_compare_means() +   # Add global p-value
    scale_fill_manual(values = c('#ffca7b','#820a86', "grey")) +
    scale_color_manual(values = c('#ffca7b','#820a86', "grey")) +
    labs(title="Effect of TEs over gene-TEs", x="", y ="Number of teM epiallele per gene")+
    theme_bw() +
    facet_grid(.~DMRoverTEs, scale="free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    theme(panel.grid.major = element_line(color = "gray90", linetype = "dashed"))
  ggsave(paste0(outPlot,"03.03.gene-TE_TEeffect_Epiallele_frequency_boxplot.pdf"),  width = 20, height = 20, units = "cm")  
  
  # Wilcoxon rank-sum test
  # remove one epiallele
  wilcoxonRes<-dmrOverTeBins %>%
    filter(sampleEpiallele=="gbM") %>%
    mutate(DMRoverTEs = factor(DMRoverTEs)) %>%
    group_by(frequency_bin) %>%
    group_modify(~ {
      test_result <- wilcox.test(count ~ DMRoverTEs, data = .x, p.adjust.method = "bonferroni")
      tidy(test_result)
    })
}
########################################################################################
##  D.2 | TE effect over genes: gene case, distance to TEs
######################################################################################## 
{
  #genes<-unique(subset(allData, geneType=="gene")$geneName)
  allData[, relativeDis := ifelse(absDistanceTE == 0, 0, 
                                  ifelse(absDistanceTE <= 100, 100, 
                                         ifelse(absDistanceTE > 100 & absDistanceTE <= 250, 250,
                                                ifelse(absDistanceTE > 250 & absDistanceTE <= 500, 500,
                                                       ifelse(absDistanceTE > 500 & absDistanceTE <= 1000,1000,1001)))))]
  allData[, relativeDistoGene := ifelse(TEposition == 'upstream', relativeDis*-1,relativeDis)]
  relativeDis<-unique(allData$relativeDis)
  distance2te<-data.table()
  for(d in relativeDis){
    tmp <- allData[relativeDis == d ]
    summ<-tmp %>%
      group_by(geneName, geneType, geneEpiallele, annotation,TEposition,absDistanceTE,relativeDis,relativeDistoGene,sampleEpiallele) %>%
      dplyr::summarise(count = n(),frequency = n()/186, .groups = 'drop') %>%
      ungroup()
    distance2te<- rbindlist(list(distance2te, summ), fill = TRUE)
  }
  distance2te
  tmp <- distance2te %>%
    filter(sampleEpiallele=="teM")
  tmp <- distance2te
  tmp$relativeDistoGene <- factor(tmp$relativeDistoGene, levels = sort(unique(tmp$relativeDistoGene)))
  tmp$relativeDis <- factor(tmp$relativeDis, levels = sort(unique(tmp$relativeDis)))
  
  ggplot(tmp, aes(x = relativeDistoGene , y = frequency, fill = sampleEpiallele)) +
    geom_boxplot(width=0.5, na.rm = T,lwd = 0.3,alpha=0.7,outlier.shape = NA)+
    #stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred") + 
    labs(x = "TE distance to the closest PE gene", y = "Sample frequency")#+
  #scale_fill_manual(values =  '#ffca7b') +
  #scale_color_manual(values = '#ffca7b') 
  ggsave(paste0(outPlot,"03.04.gene_TEeffect_epiallele-frequency_close-far_boxplot.pdf"),  width = 20, height = 20, units = "cm")
  
  # Perform the Wilcoxon rank-sum test
  toCompare<-c("100","1001")
  compare100to1001 <- tmp %>%
    filter(relativeDis %in% toCompare)
  wilcox.test(
    frequency ~ relativeDis,
    data = compare100to1001
  )
  # Wilcoxon rank sum test with continuity correction
  toCompare<-c("-100","-1001")
  compare100to1001 <- tmp %>%
    filter(relativeDistoGene %in% toCompare)
  wilcox.test(
    frequency ~ relativeDis,
    data = compare100to1001
  )
  
  ggplot(tmp, aes(x = absDistanceTE, y = frequency)) +
    geom_point() +
    labs(x = "teDistance bin", y = "Frequency") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle("Scatter Plot of Frequency vs. teDistance bin")+
    facet_grid(.~TEposition)
  ggsave(paste0(outPlot,"03.05.gene_TEeffect_epiallele-frequency_close-far_dot.pdf"),  width = 20, height = 20, units = "cm")
}
