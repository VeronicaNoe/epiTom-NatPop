########################################################################################
##  Set directories and load files
######################################################################################## 
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(broom)
  library(ggplot2)
  
  library(gplots)
  library(ggrepel)
  library(dendextend)
  library(sva) # combat
  library(DESeq2)
  library(ggpubr)
  library(plyr)
  
  library(purrr)
  library(FSA)
  
  library(FactoMineR)
  library(factoextra)
  
  library(data.table)
  library("RColorBrewer")
  
  library(viridis)
  
})
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig4/results")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/plots/"

{
geneExp<-data.table::fread("/mnt/disk2/vibanez/07_rnaseq-processing/07.3_counts/leaf-transcriptome_deseq-combat.tsv", sep = '\t', 
                           data.table = T, fill = TRUE, check.names=FALSE,
                           na.string=c("NA"), nThread = 10)
geneExp <- geneExp %>%
    pivot_longer(
      cols = -Geneid,     
      names_to = "Sample",  
      values_to = "Expression"    
    )
  samplesWexp<-unique(geneExp$Sample)
  geneExp$expressionLog2<-log2(geneExp$Expression+1)

allData<-data.table::fread('01.06_epiallele_sample.tsv', sep = '\t', 
                           data.table = T, fill = TRUE, check.names=FALSE,
                           na.string=c("NA"), nThread = 20)

bySamplesEpiallele <- allData %>%
  filter(accessions %in% samplesWexp)
gc()

bySamplesEpiallele <- bySamplesEpiallele  %>%
  left_join(geneExp, by = c("geneName" = "Geneid", "accessions"= "Sample"))
rm(allData)  
rm(geneExp)  
gc()
data.table::fwrite(bySamplesEpiallele, file = paste0(outDir, "04.01.epiallele-geneExpression.tsv"))
}
########################################################################################
##  A |  Expression of genes affected by meth state
######################################################################################## 
{
bySamplesEpiallele<-data.table::fread(paste0(outDir,"04.01.epiallele-geneExpression.tsv")
                                        , sep = ',', 
                             data.table = T, fill = TRUE, check.names=FALSE,
                             na.string=c("NA"), nThread = 10)
  
length(unique(bySamplesEpiallele$geneName))
# get genes with 2 or 3 epialleles
informativeGenes <- bySamplesEpiallele %>%
  select(geneName, geneEpiallele,annotation, accessions, sampleEpiallele, expressionLog2) %>%
  group_by(geneName,accessions) %>%
  dplyr::summarise(geneEpiallele = unique(geneEpiallele),
                   ann = paste(unique(annotation), collapse = ","),
                   sampleEpiallele =paste(unique(sampleEpiallele), collapse = ","),
                   expressionLog2 = mean(expressionLog2)
                   ,.groups = 'drop') %>%
  group_by(geneName, ann) %>%
  dplyr::summarise(
    meanExpression = mean(expressionLog2, na.rm = TRUE),
    nSampleEpi = paste(unique(sampleEpiallele), collapse = ","),
    minEpi = length(unique(sampleEpiallele)),
    nSamples = n(),.groups = 'drop') %>%
  filter( minEpi >1)

#  keep genes with more than 3 samples per epiallele per gene
informativeData<-bySamplesEpiallele %>%
  filter(geneName %in% informativeGenes$geneName) %>%
  group_by(geneName,sampleEpiallele) %>%
  dplyr::summarise(nSamples=n()) %>%
  group_by(geneName) %>%
  dplyr::summarise(minSample=min(nSamples))%>%
  filter(minSample>3)

# filter data
dt<-bySamplesEpiallele %>%
  filter(geneName %in% informativeData$geneName) %>%
  select(geneName,annotation,accessions,sampleEpiallele,expressionLog2)
length(unique(dt$geneName))
kruskal_results <- dt %>%
  group_by(geneName) %>%
  do(tidy(kruskal.test(expressionLog2 ~ sampleEpiallele, data = .))) %>%
  ungroup()
length(unique(subset(kruskal_results,p.value<0.05)$geneName))


genes2plot<-dt %>%
  select(geneName,annotation,accessions,expressionLog2,sampleEpiallele)

genes2plot <- genes2plot %>%
  left_join(kruskal_results %>% select(geneName, p.value), by = "geneName")

samplEpi_gene<-genes2plot %>%
  group_by(geneName) %>%
  dplyr::summarise(sEpi = paste(sort(unique(sampleEpiallele)), collapse = ","),
                   nEpi= ifelse(sEpi=="gbM,UM" |sEpi=="teM,UM",'biEpiAllele','polyEpiAllele' ))

genes2plot <- genes2plot %>%
  inner_join(samplEpi_gene, by = "geneName") %>%
  filter(sEpi!='gbM,teM')

# change direction
gbM_bi <- genes2plot %>%
  filter(sEpi == "gbM,UM")%>%
  group_by(geneName, sampleEpiallele,p.value) %>%
  dplyr::summarise(meanExp = mean(expressionLog2)) %>%
  pivot_wider(names_from = sampleEpiallele, values_from = meanExp) %>%
  group_by(geneName,UM,gbM,p.value) %>%
  dplyr::summarise(logFC = gbM-UM,
    change = ifelse(gbM-UM<0,'negative','positive'))

merged_gbM <- gbM_bi %>%
  mutate(logP = -log10(p.value)) %>%
  mutate(significance = case_when(
    logFC >= 1 & logP > -log10(0.05) ~ "upregulated",
    logFC <= -1  & logP > -log10(0.05) ~ "downregulated",
    (logFC < 1 | logFC > -1 )& logP > -log10(0.05) ~ "moderated",
    TRUE ~ "neutral"
  ))
data.table::fwrite(merged_gbM, file = paste0(outDir, "04.02.gbM-geneExpression.tsv"))
nDEG<-nrow(subset(merged_gbM, significance=="upregulated" | significance=="downregulated"))
# Plot the volcano plot
ggplot(merged_gbM, aes(x = logFC, y = logP, color=significance,fill=significance)) +
  geom_point(size=3) +
  labs(title = paste0(nDEG," DEG UM vs gbM epialleles"),
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")+
  scale_color_manual(values = c("upregulated" = "#01cda560", "downregulated" = "#bb0d0085", 
                                "moderated" ="grey","neutral"="grey"), 
                     name = "Regulation")+
  scale_fill_manual(values = c("upregulated" = "#01cda560", "downregulated" = "#bb0d0085", 
                              "moderated" ="grey","neutral"="grey"), 
                   name = "Regulation")
ggsave(paste0(outPlot,"04.01.geneExpression-gbM-Epiallele_volcano.pdf"),  width = 20, height = 20, units = "cm")
merged_gbM %>%
  group_by(significance) %>%
  dplyr::summarise(n())
# A tibble: 4 × 2
# significance  `n()`
# <chr>         <int>
# 1 downregulated    65
# 2 moderated       685
# 3 neutral       11242
# 4 upregulated      84
###
teM_bi <- genes2plot %>%
  filter(sEpi == "teM,UM")%>%
  group_by(geneName, sampleEpiallele,p.value) %>%
  dplyr::summarise(meanExp = mean(expressionLog2)) %>%
  pivot_wider(names_from = sampleEpiallele, values_from = meanExp) %>%
  group_by(geneName,UM,teM,p.value) %>%
  dplyr::summarise(logFC = teM-UM,
                   change = ifelse(teM-UM<0,'negative','positive'))

merged_teM <- teM_bi %>%
  mutate(logP = -log10(p.value)) %>%
  mutate(significance = case_when(
    logFC >= 1 & logP > -log10(0.05) ~ "upregulated",
    logFC <= -1  & logP > -log10(0.05) ~ "downregulated",
    (logFC < 1 | logFC > -1 )& logP > -log10(0.05) ~ "moderated",
    TRUE ~ "neutral"
  ))
data.table::fwrite(merged_teM, file = paste0(outDir, "04.3.teM-geneExpression.tsv"))
nDEG<-nrow(subset(merged_teM, significance=="upregulated" | significance=="downregulated"))
# Plot the volcano plot
ggplot(merged_teM, aes(x = logFC, y = logP, color=significance)) +
  geom_point(size=3) +
  labs(title = paste0(nDEG," UM vs teM epialleles"),
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")+
  scale_color_manual(values = c("upregulated" = "#01cda560", "downregulated" = "#bb0d0085", 
                                "moderated" ="grey","neutral"="grey"), 
                     name = "Regulation")+
  scale_fill_manual(values = c("upregulated" = "#01cda560", "downregulated" = "#bb0d0085", 
                               "moderated" ="grey","neutral"="grey"), 
                    name = "Regulation")+
  theme_minimal()
ggsave(paste0(outPlot,"04.02.geneExpression-teM-Epiallele_volcano.pdf"),  width = 20, height = 20, units = "cm")
merged_teM %>%
  group_by(significance) %>%
  dplyr::summarise(n())
# # A tibble: 4 × 2
# significance  `n()`
# <chr>         <int>
#   1 downregulated   389
# 2 moderated       403
# 3 neutral        4479
# 4 upregulated      47
#
# Create the plot with counts and p-values
merged_teM %>%
  filter(significance=="downregulated" & logP>20)


##pe
pe <- genes2plot %>%
  filter(sEpi == "gbM,teM,UM")%>%
  group_by(geneName, sampleEpiallele,p.value) %>%
  dplyr::summarise(meanExp = mean(expressionLog2)) %>%
  pivot_wider(names_from = sampleEpiallele, values_from = meanExp) 
  
significant_pe <- pe %>%
  filter(p.value < 0.05)%>%
  mutate(
    UM_gbM_change = gbM - UM,
    UM_teM_change = teM - UM,
    significant_direction = case_when(
      UM_gbM_change > 0 & UM_teM_change < 0 ~ "Up-gbM_Down-teM",
      UM_gbM_change < 0 & UM_teM_change > 0 ~ "Down-gbM_Up-teM",
      UM_gbM_change > 0 & UM_teM_change > 0 ~ "Up",
      UM_gbM_change < 0 & UM_teM_change < 0 ~ "Down",
      TRUE ~ "non"
    )
  )
data.table::fwrite(significant_pe, file = paste0(outDir, "04.04.PE-geneExpression.tsv"))
long_significant_genes <- significant_pe %>%
  filter((UM_gbM_change< -1 | UM_gbM_change > 1) | (UM_teM_change < -1 |UM_teM_change >1))%>%
  pivot_longer(cols = c("UM", "gbM", "teM"), 
               names_to = "sampleEpiallele", 
               values_to = "expressionLog2")

nDEG<-significant_pe %>%
             filter((UM_gbM_change< -1 | UM_gbM_change > 1) | (UM_teM_change < -1 |UM_teM_change >1))
# Plot the data
ggplot(long_significant_genes, aes(x = sampleEpiallele, y = expressionLog2, 
                                   color = sampleEpiallele, fill = sampleEpiallele)) +
  geom_boxplot(width = 0.5, na.rm = TRUE, alpha = 0.6, lwd = 0.3, outlier.shape = NA) +
  geom_line(aes(group = geneName),color="grey90", lwd = 0.5, alpha = 0.8, na.rm = TRUE) +  # Lines connecting points by gene
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", 
               position = position_dodge(width = 0.75)) +
  geom_point(size = 1.5, alpha = 0.7, na.rm = TRUE) +  # Add jittered points
  scale_fill_manual(values = c('#ffca7b', '#820a86', "grey")) +
  scale_color_manual(values = c('#ffca7b','#820a86', "grey")) +
  labs(title = paste0(nrow(nDEG)," DEG for PE genes"),
       x = "Epiallele Type", y = "Expression Log2") +
  theme_bw() +
  facet_wrap(~ significant_direction) +
  theme(panel.grid.major = element_line(color = "gray90", linetype = "dashed"))
ggsave(paste0(outPlot,"04.03.geneExpression-PolyEpiallele_boxplot.pdf"),  width = 20, height = 20, units = "cm")

nDEG %>%
  group_by(significant_direction) %>%
  dplyr::summarise(count = n())
}
# A tibble: 3 × 2
# significant_direction count
# <chr>                 <int>
# 1 Down                     28
# 2 Up                        2
# 3 Up-gbM_Down-teM          11
long_significant_genes %>%
  filter(significant_direction=="Up-gbM_Down-teM") %>%
  select(geneName)
########################################################################################
##  B | heatmap per gene - sample
######################################################################################## 
{## get wider table
  
  gene2heat<-genes2plot %>%
    group_by(geneName,sampleEpiallele, p.value, sEpi,nEpi) %>%
    dplyr::summarise(meanExp = mean(expressionLog2))
  
### get epi information
epiTab<-data.frame(epialleleType = c("gbM", "teM", "UM"),
                   colors= c("#ffca7b", "#820a86", "grey"))

### get color for expression
colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
my_palette <- colorRampPalette(c("#00136f", "#9e0039", "#136f00"))(n = 299)
nEpiAlleles<-unique(gene2heat$nEpi)
#nEpiAlleles<-"polyEpiAllele"

for(nE in nEpiAlleles){
    tmp<-subset(gene2heat, nEpi==nE)
    epiTypes<-unique(tmp$sEpi)
    for(tE in epiTypes){
      toHeat<-as.data.frame(subset(tmp, sEpi==tE))
      toHeat<-toHeat %>%
        filter(sEpi==tE) %>%
        pivot_wider(names_from = sampleEpiallele, values_from = meanExp)
      toHeat<-as.data.frame(toHeat)
      toHeat$p.value<-NULL
      rownames(toHeat)<-toHeat$geneName
      toHeat$sEpi<-NULL
      toHeat$nEpi<-NULL
      toHeat$geneName<-NULL
      
      
      toGroup<-epiTab[epiTab$epialleleType %in% colnames(toHeat),]
      toGroup<-as.data.frame(toGroup)
      rownames(toGroup)<-toGroup$epialleleType
      toGroup<-toGroup[colnames(toHeat),]
      
      gene_dist <- dist(toHeat)  # Calculate distance matrix
      gene_hclust <- hclust(gene_dist)  # Perform hierarchical clustering
      gene_dendrogram <- as.dendrogram(gene_hclust)  # Convert to dendrogram
      
      cut_height <- 20  # Adjust this value as needed
      
      clusters <- cutree(gene_hclust, h = cut_height)
      gene_clusters <- split(names(clusters), clusters)
      gene_clusters_df <- data.frame(
        geneName = unlist(gene_clusters, use.names = FALSE),
        cluster = rep(names(gene_clusters), times = sapply(gene_clusters, length)))
      
      #data.table::fwrite(gene_clusters_df, file = paste0(outDir, "05.1.a.clusters_", g,"_", a,".csv"))
      gene_dendrogram <- color_branches(gene_dendrogram, h = cut_height)  #
      #plot(gene_dendrogram)
      gc()
      ## do the heatmap plot
      pdf(paste0(outPlot, "03.4.", nE, "_",tE,"_geneExpression_heatmap.pdf"), width = 20, height = 20)
      heatmap.2(as.matrix(toHeat), scale = "none",
                labCol = toGroup[,"epialleleType"],
                #labRow = geneInfo2plot[,"geneEpiallele"],
                #RowSideColors = geneInfo2plot[,"epialleleCol"],
                ColSideColors = as.character(toGroup[,"colors"]),
                col = my_palette,
                trace = "none", density.info = "none",
                Rowv = gene_dendrogram)  # Add clustering dendrogram
      dev.off()
      gc()
    }
  }
}
########################################################################################
##  C | Pe gene expression f(cambio de frecuencia)
######################################################################################## 
geneList<-data.table::fread(paste0(outDir,"02.02.accGroups_teM-frequency_PIM-SLC.tsv"), 
                            sep = '\t', header=T,
                              data.table = T, fill = TRUE, check.names=FALSE,
                              na.string=c("NA"), nThread = 10)

geneList<-geneList %>%
  filter(p.adjusted_fdr<=0.05)

length(unique(geneList$geneName))
merged_teM %>%
  filter(geneName %in% geneList$geneName ) %>%
  group_by(significance)%>%
  dplyr::summarise(nCases= n())

# A tibble: 4 × 2
# significance  nCases
# <chr>          <int>
# 1 downregulated     55
# 2 moderated        104
# 3 neutral          850
# 4 upregulated       12
c# A tibble: 4 × 2
# significant_direction nCases
# <chr>                  <int>
# 1 Down                      25
# 2 Down-gbM_Up-teM            2
# 3 Up                         5
# 4 Up-gbM_Down-teM            8
DEG<-significant_pe %>%
  filter(geneName %in% geneList$geneName ) %>%
  group_by(significant_direction)
DEG<-unique(DEG$geneName)

tmp<-geneList %>%
  filter(geneName %in% DEG) 
PimSlc <- tmp %>%
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
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
   geom_text_repel(data = PimSlc %>% filter(sig_levels == 0.00001), 
                      aes(label = geneName), size = 2, max.overlaps = 50) +  # Add gene labels for significant points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Perfect correlation line
  scale_color_identity() +
  scale_fill_identity() +
  #facet_grid(geneType ~ .)+
  labs(x = "PIM", y = "SLC") +
  theme_bw()
ggsave(paste0(outPlot,"04.04.geneExpression-PimSlc.pdf"),  width = 20, height = 20, units = "cm")
