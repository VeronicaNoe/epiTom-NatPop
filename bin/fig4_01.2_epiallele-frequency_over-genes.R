suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(ggpubr)
  library(data.table)
  library(tidyr)
  library(purrr)
})
setwd("10_data-analysis/Fig4/aa_get-epialleles-over-genes/bd_merged-gene-epialleles")
outDir<-"10_data-analysis/Fig4/results/"
outPlot<-"10_data-analysis/Fig4/plots/"


allData<-data.table::fread(paste0(outDir,'01.06_epiallele_sample.tsv'), sep = '\t',
                           data.table = T, fill = TRUE, check.names=FALSE,
                           na.string=c("NA"), nThread = 10)
########################################################################################
##  A |  epiallelic frequency over genes
######################################################################################## 
## Get the histogram of frequencies
{
  allele_frequency <- allData %>%
    group_by(geneName,geneType,geneEpiallele,annotation, sampleEpiallele,teMnoverTEs) %>%
    dplyr::summarise(count = n()/186) %>%
    mutate(freqClass = ifelse(count>0.9,'a.highFreq',
                              ifelse(count>0.05 & count<=0.9,'b.intermediateFreq',
                                     'c.lowFreq')))
  
  data.table::fwrite(allele_frequency, file=paste0(outDir,"02.01.epiallele-frequency.tsv"),
                     quote=F,row.names=T,col.names = T,sep="\t")
  
  ggplot(allele_frequency, aes(x = count, fill = sampleEpiallele, color=sampleEpiallele)) +
    geom_histogram(aes(y = ..density..), bins = 100, position = "identity", alpha = 0.5) +  # Density on y-axis
    scale_fill_manual(values = c('#ffca7b', '#820a86', "grey")) +
    scale_color_manual(values = c('#ffca7b', '#820a86', "grey")) +
    labs(title = "Density Histogram of Allele Frequencies",
         x = "Allele Frequency",
         y = "Density") +  # Label change to reflect density
    scale_x_continuous(breaks = seq(0, max(allele_frequency$count), by = 0.05)) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(paste0(outPlot,"02.01.nEpiallele_frequency_histogram.pdf"),  width = 20, height = 20, units = "cm")  
}
## Get the boxplots of frequencies
{
  toPlot<- allData %>%
    group_by(geneName, geneType, sampleEpiallele, annotation) %>%
    dplyr::summarise(count = n()) %>%
    ungroup()
  
  ggplot(toPlot, aes(x = sampleEpiallele, y = count, color=sampleEpiallele, fill=sampleEpiallele)) +
    geom_boxplot(width=0.5, na.rm = T,alpha=0.6,lwd = 0.3,outlier.shape = NA)+
    stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", position = position_dodge(width = 0.75)) +
    geom_violin(width = 0.7, alpha = 0.4, lwd = 0.3, position = position_dodge(width = 0.75)) +
    stat_compare_means() +   # Add global p-value
    scale_fill_manual(values = c('#ffca7b', '#820a86', "grey")) +
    scale_color_manual(values = c('#ffca7b', '#820a86', "grey")) +
    labs(title="", x="", y ="Number of epiallele per gene")+
    theme_bw() +
    facet_grid(. ~ geneType, scales = 'free') +
    theme(panel.grid.major = element_line(color = "gray90", linetype = "dashed"))
  ggsave(paste0(outPlot,"02.02.nEpiallele_frequency_boxplot.pdf"),  width = 20, height = 20, units = "cm")
}
########################################################################################
##  B | TE frequency among tomato groups
######################################################################################## 
{
  setDT(allData)
# get plots of global frequency per groups
{
  groupsFreq <- allData %>%
    group_by(geneName, geneEpiallele, geneType, annotation, sampleEpiallele,accGroup) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    group_by(geneName, geneEpiallele, geneType, annotation, accGroup) %>%
    dplyr::mutate(total_count_accGroup = sum(count)) %>%  ## don't forget dplyr:: it's not the same without that
    ungroup() %>%
    mutate(frequency = count / total_count_accGroup) 
  
  groupsFreq$accGroup <- factor(groupsFreq$accGroup, levels = c("WILD", "PIM", "SLC", "SLL"))
  ggplot(groupsFreq, aes(x = frequency, fill = accGroup, color=accGroup, alpha=0.4)) +
    geom_histogram(position = "identity", breaks = seq(0, max(groupsFreq$frequency), by = 0.04)) +
    scale_fill_manual(values = c("#000000ff", "#00a681ff", "#ed83b5ff", "#00a4deff")) +
    scale_color_manual(values = c("#000000ff", "#00a681ff", "#ed83b5ff", "#00a4deff")) +
    scale_x_continuous(breaks = seq(0, max(groupsFreq$frequency), by = 0.1)) +  
    labs(title = "Histogram of Allele Frequencies",
         x = "Allele Frequency",
         y = "Number of gene per epiallele") +
    facet_grid(sampleEpiallele~., scale="free") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  ggsave(paste0(outPlot,"02.03.accGroups_Epiallele_frequency_histogram.pdf"),  width = 20, height = 20, units = "cm")  
  
  groupsFreqBins <- groupsFreq %>%
    mutate(frequency_bin = cut(frequency, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE))
  # Count the occurrences within each bin
  frequency_counts <- groupsFreqBins %>%
    #filter(sampleEpiallele!="gbM") %>%
    group_by(geneType, geneEpiallele, sampleEpiallele,annotation, accGroup, frequency_bin) %>%
    dplyr::summarize(count = n())
  frequency_counts$accGroup <- factor(frequency_counts$accGroup, levels = c("WILD", "PIM", "SLC", "SLL"))
  ## Get the boxplots of frequencies
    ggplot(frequency_counts, aes(x = frequency_bin, y = count, color=accGroup, 
                                 fill=accGroup)) +
      geom_boxplot(width=0.5, na.rm = F,alpha=0.6,lwd = 0.3,outlier.shape = NA)+
      stat_summary(fun = mean, geom = "point", shape = 16, size = 2, 
                   color = "darkred", position = position_dodge(width = 0.5)) + 
      #stat_compare_means() +   # Add global p-value
      scale_fill_manual(values = c("#000000ff", "#00a681ff", "#ed83b5ff", "#00a4deff")) +
      scale_color_manual(values = c("#000000ff", "#00a681ff", "#ed83b5ff", "#00a4deff")) +
      labs(title="Epiallele frequency per group", x="", y ="")+
      theme_bw() +
      facet_grid(sampleEpiallele~., scale="free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
      theme(panel.grid.major = element_line(color = "gray90", linetype = "dashed"))
    ggsave(paste0(outPlot,"02.04.accGroups_Epiallele_frequency_boxplot.pdf"),  width = 20, height = 20, units = "cm")  
}
# compare frequency per group
  groupsFreq
  teMfreq <- groupsFreq %>%
    filter(sampleEpiallele=="teM") %>%
    pivot_wider(
      id_cols = c(geneName,geneType,geneEpiallele,annotation),
      names_from = accGroup,
      values_from = frequency,
      values_fill = list(frequency = 0))
  ##
  total_samples_PIM <- 38 
  total_samples_SLC <- 75 
  total_samples_SLL <- 50 
  ## PIM-SLC
  {
  PimSlc <- teMfreq %>%
    mutate(PIM_count = round(PIM * total_samples_PIM),
      SLC_count = round(SLC * total_samples_SLC))
  # func
  test_PIM_SLC <- function(pim_count, slc_count, total_pim, total_slc) {
    if (pim_count == 0 && slc_count == 0) {
      return(NA)
    }
    contingency_table <- matrix(c(pim_count, total_pim - pim_count, slc_count, total_slc - slc_count), nrow = 2)
    if (min(contingency_table) < 5) {
      test_result <- tryCatch({
        fisher.test(contingency_table)$p.value
      }, error = function(e) NA)
    } else {
      test_result <- chisq.test(contingency_table, simulate.p.value = TRUE)$p.value
    }
    return(test_result)
  }
  # Apply the test for each gene using pmap
  PimSlc <- PimSlc %>%
    mutate(p_value = pmap_dbl(list(PIM_count, SLC_count, total_samples_PIM, total_samples_SLC), test_PIM_SLC))%>%
    mutate(p.adjusted_bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
    mutate(p.adjusted_fdr = p.adjust(p_value, method = "fdr"))
  
  data.table::fwrite(PimSlc, file=paste0(outDir,"02.02.accGroups_teM-frequency_PIM-SLC.tsv"),
                     quote=F,row.names=T,col.names = T,sep="\t")
  
  PimSlc <- PimSlc %>%
    mutate(sig_levels = ifelse(p.adjusted_bonferroni <= 0.05 & p.adjusted_bonferroni > 0.001, 0.05,
                               ifelse(p.adjusted_bonferroni <= 0.001 & p.adjusted_bonferroni > 0.00001, 0.001,
                                      ifelse(p.adjusted_bonferroni <=0.00001, 0.00001,0))))
  # for plot
  PimSlc <- PimSlc %>%
    mutate(sig_color = case_when(
      sig_levels == 0.05 ~  "#C2B96B",
      sig_levels == 0.001 ~ "#e69505",
      sig_levels == 0.00001 ~"#ed1002",
      TRUE ~ "grey"
    ))
  
  set.seed(1234)
  # Create the plot with the new color column and add a legend
  ggplot(PimSlc, aes(x = PIM, y = SLC, color = sig_color, fill = sig_color)) +
    geom_jitter(size = 3, shape = 21, width = 0.01, height = 0.01) +  
    #geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
    # geom_text_repel(data = PimSlc %>% filter(sig_levels == 0.00001), 
    #                    aes(label = geneName), size = 2, max.overlaps = 50) +  # Add gene labels for significant points
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Perfect correlation line
    scale_color_identity() +
    scale_fill_identity() +
    #facet_grid(geneType ~ .)+
    labs(x = "PIM", y = "SLC") +
    theme_bw()
  ggsave(paste0(outPlot,"02.05.teM-frequency_PIM-SLC_bonferroni.pdf"),  width = 20, height = 20, units = "cm")
  
  # fdr
  PimSlc <- PimSlc %>%
    mutate(sig_levels = ifelse(p.adjusted_fdr <= 0.05 & p.adjusted_fdr > 0.001, 0.05,
                               ifelse(p.adjusted_fdr <= 0.001 & p.adjusted_fdr > 0.00001, 0.001,
                                      ifelse(p.adjusted_fdr <=0.00001, 0.00001,0))))
  # for plot
  PimSlc <- PimSlc %>%
    mutate(sig_color = case_when(
      sig_levels == 0.05 ~  "#C2B96B",
      sig_levels == 0.001 ~ "#e69505",
      sig_levels == 0.00001 ~"#ed1002",
      TRUE ~ "grey"
    ))
  ggplot(PimSlc, aes(x = PIM, y = SLC, color = sig_color, fill = sig_color)) +
    geom_jitter(size = 3, shape = 21, width = 0.01, height = 0.01) +  
    #geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
    # geom_text_repel(data = PimSlc %>% filter(sig_levels == 0.00001), 
    #                    aes(label = geneName), size = 2, max.overlaps = 50) +  # Add gene labels for significant points
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Perfect correlation line
    scale_color_identity() +
    scale_fill_identity() +
    #facet_grid(geneType ~ .)+
    labs(x = "PIM", y = "SLC") +
    theme_bw()
  ggsave(paste0(outPlot,"02.06.teM-frequency_PIM-SLC_fdr.pdf"),  width = 20, height = 20, units = "cm")
  
  }
  ## SLC-SLL
  {
  SlcSll <- teMfreq %>%
    mutate(SLC_count = round(SLC * total_samples_SLC),
           SLL_count = round(SLL * total_samples_SLL))
  # func
  test_SLC_SLL <- function(slc_count,sll_count, total_slc,total_sll) {
    if (sll_count == 0 && slc_count == 0) {
      return(NA)
    }
    contingency_table <- matrix(c(slc_count, total_slc - slc_count,
                                  sll_count, total_sll - sll_count), nrow = 2)
    if (min(contingency_table) < 5) {
      test_result <- tryCatch({
        fisher.test(contingency_table)$p.value
      }, error = function(e) NA)
    } else {
      test_result <- chisq.test(contingency_table, simulate.p.value = TRUE)$p.value
    }
    return(test_result)
  }
  # Apply the test for each gene using pmap
  SlcSll <- SlcSll %>%
    mutate(p_value = pmap_dbl(list(SLC_count, SLL_count, total_samples_SLC, 
                                   total_samples_SLL), test_SLC_SLL))%>%
    mutate(p.adjusted_bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
    mutate(p.adjusted_fdr = p.adjust(p_value, method = "fdr"))
  
  data.table::fwrite(SlcSll, file=paste0(outDir,"02.03.accGroups_teM-frequency_SLC-SLL.tsv"),
                     quote=F,row.names=T,col.names = T,sep="\t")
  
  SlcSll <- SlcSll %>%
    mutate(sig_levels = ifelse(p.adjusted_bonferroni <= 0.05 & p.adjusted_bonferroni > 0.001, 0.05,
                               ifelse(p.adjusted_bonferroni <= 0.001 & p.adjusted_bonferroni > 0.00001, 0.001,
                                      ifelse(p.adjusted_bonferroni <=0.00001, 0.00001,0))))
  # for plot
  SlcSll <- SlcSll %>%
    mutate(sig_color = case_when(
      sig_levels == 0.05 ~  "#C2B96B",
      sig_levels == 0.001 ~ "#e69505",
      sig_levels == 0.00001 ~"#ed1002",
      TRUE ~ "grey"
    ))
  
  set.seed(1234)
  # Create the plot with the new color column and add a legend
  ggplot(SlcSll, aes(x = SLC, y = SLL, color = sig_color, fill = sig_color)) +
    geom_jitter(size = 3, shape = 21, width = 0.01, height = 0.01) +  
    #geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
    # geom_text_repel(data = PimSlc %>% filter(sig_levels == 0.00001), 
    #                    aes(label = geneName), size = 2, max.overlaps = 50) +  # Add gene labels for significant points
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Perfect correlation line
    scale_color_identity() +
    scale_fill_identity() +
    #facet_grid(geneType ~ .)+
    labs(x = "SLC", y = "SLL") +
    theme_bw()
  ggsave(paste0(outPlot,"02.07.teM-frequency_SLC-SLL_bonferroni.pdf"),  width = 20, height = 20, units = "cm")
  
  # fdr
  SlcSll <- SlcSll %>%
    mutate(sig_levels = ifelse(p.adjusted_fdr <= 0.05 & p.adjusted_fdr > 0.001, 0.05,
                               ifelse(p.adjusted_fdr <= 0.001 & p.adjusted_fdr > 0.00001, 0.001,
                                      ifelse(p.adjusted_fdr <=0.00001, 0.00001,0))))
  # for plot
  SlcSll <- SlcSll %>%
    mutate(sig_color = case_when(
      sig_levels == 0.05 ~  "#C2B96B",
      sig_levels == 0.001 ~ "#e69505",
      sig_levels == 0.00001 ~"#ed1002",
      TRUE ~ "grey"
    ))
  ggplot(SlcSll, aes(x = SLC, y = SLL, color = sig_color, fill = sig_color)) +
    geom_jitter(size = 3, shape = 21, width = 0.01, height = 0.01) +  
    #geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
    # geom_text_repel(data = PimSlc %>% filter(sig_levels == 0.00001), 
    #                    aes(label = geneName), size = 2, max.overlaps = 50) +  # Add gene labels for significant points
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Perfect correlation line
    scale_color_identity() +
    scale_fill_identity() +
    #facet_grid(geneType ~ .)+
    labs(x = "SLC", y = "SLL") +
    theme_bw()
  ggsave(paste0(outPlot,"02.08.teM-frequency_SLC-SLL_fdr.pdf"),  width = 20, height = 20, units = "cm")
  
  ### to check
  # tmp<-PimSlc %>% filter(sig_levels==0.001 & SLC>=0.7)
  # tmp$geneName
  # #high SLC
  # # Solyc10g076960.1
  # subset(tmp, geneName =="Solyc06g054370.1")
  # 
  # subset(allData,geneName =="Solyc06g054370.1" & accessions=="TS-95")
  # epialleleState<- data.table::fread(paste0(outDir,"allSamples_epiAlleles.tsv"), sep = '\t',
  #                                    data.table = T, fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
  # subset(epialleleState,geneName =="Solyc02g068820.1" &
  #          accessions=="TS-95" )
}
  ### merge comparisons
  {
  SlcSll_filtered <- SlcSll %>%
    select(geneName,geneType, geneEpiallele, annotation, SLL,PIM,SLC,WILD,p.adjusted_bonferroni, p.adjusted_fdr)
  PimSlc_filtered <- PimSlc %>%
    select(geneName,geneType, geneEpiallele, annotation, SLL,PIM,SLC,WILD,p.adjusted_bonferroni, p.adjusted_fdr) 
  merged_data <-PimSlc_filtered  %>%
    full_join(SlcSll_filtered, by = c("geneName","geneType", "geneEpiallele", "annotation", "SLL","PIM","SLC","WILD"), suffix = c("", ".SlcSll"))
  
  merged_data <- merged_data %>%
    mutate(
      comparison = case_when(
        p.adjusted_bonferroni < 0.05 & p.adjusted_bonferroni.SlcSll < 0.05 ~ paste0("doble:", PIM < SLC, ":", SLC < SLL),
        p.adjusted_bonferroni < 0.05 & p.adjusted_bonferroni.SlcSll > 0.05  ~ paste0("single:upSLC:", PIM < SLC),
        p.adjusted_bonferroni > 0.05 & p.adjusted_bonferroni.SlcSll < 0.05 ~ paste0("single:upSLL:", SLC < SLL),
        TRUE ~ "not significant"
      )
    )
  merged_data %>%
    dplyr::count(comparison)
  
  comp<-grep("not significant", unique(merged_data$comparison), invert = T, value=T)
  for (c in comp){
    name2save<-gsub(':',"_", c)
    tmp<-subset(merged_data, comparison==c)
    data.table::fwrite(as.data.frame(tmp$geneName), file=paste0(outDir,"04.3.",name2save,"_geneList.tsv"),
                       quote=F,row.names=F,col.names = F,sep="\t")
  }
}
