suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  #library(plyr)
  library(dplyr)
  #library(viridis)
  #library(scales)
})

###     ACHTUNG
#### AFTER CHECK FILES, CHANGE wd

setwd("/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/be_merged_correlations/old_results")
#/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/be_merged_correlations
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/plots/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/results/"
# Set options to display numbers without scientific notation
options(scipen = 999)
################################################################################
# Create a function to calculate numeric positions
calculate_numeric_position <- function(x) {
  if (grepl(",", x)) {
    # If the string contains a comma, calculate the mean of numeric parts
    parts <- as.numeric(gsub("[^0-9]+", "", strsplit(x, ",")[[1]]))
    return(mean(parts))
  } else {
    # If there's no comma, directly return the numeric part
    return(as.numeric(gsub("[^0-9]+", "", x)))
  }
}
calculate_gene_position <- function(genePart) {
  genePart <- gsub("[0-9]+", "", genePart)  # Remove numeric characters
  if ('b' %in% genePart) {
    return(20)  # Add 20 if 'b' is present
  } else if ('d' %in% genePart) {
    return(24)  # Add 24 if 'd' is present
  } else {
    return(0)   # Default value if neither 'b' nor 'd' is present
  }
}
################################################################################
mergedCorr<-data.table::fread('allCorrelation_meth-geneExpression.tsv', sep = '\t',head=F,
                                  data.table = FALSE, fill = TRUE, check.names=FALSE,
                                  na.string=c("NA"), nThread = 10)
colnames(mergedCorr)<-c("chr","start","DMR","geneID","strand","genePart","pearsonExpMeth", 
                        "pValpearsonExpMeth","spearmanExpMeth","pValSpearmanExpMeth")
head(mergedCorr)
##########################################
# rename -/+ strands
##########################################
keepA <- grep('^a[^b]*$', unique(mergedCorr$genePart), value = TRUE)
tmpA<-mergedCorr[mergedCorr$genePart %in% keepA, ]
#### this shouldn't be necessary as I already take into account the gene position to get up / down stream
# tmpApos<-subset(tmpApos, strand=="+")
# 
# tmpAneg<-mergedCorr[mergedCorr$genePart %in% keepA, ]
# tmpAneg<-subset(tmpAneg, strand=="-")
# tmpAneg$genePart<-gsub('a','d',tmpAneg$genePart)
# head(tmpAneg)

keepD <- grep('^d[^b]*$', unique(mergedCorr$genePart), value = TRUE) # keep only dmrs overlaping one category. 
tmpD<-mergedCorr[mergedCorr$genePart %in% keepD, ]
# tmpDneg<-subset(tmpDneg, strand=="-")
# tmpDpos<-mergedCorr[mergedCorr$genePart %in% keepD, ]
# tmpDpos<-subset(tmpDpos, strand=="+")
# tmpDpos$genePart<-gsub('d','a',tmpDpos$genePart)
# head(tmpDpos)

keepB <- grep('^b[^a]*$', unique(mergedCorr$genePart), value = TRUE) # keep only dmrs overlaping one category. 
keepB <- grep('^b[^d]*$', keepB, value = TRUE) # keep only dmrs overlaping one category. 
tmpB<-mergedCorr[mergedCorr$genePart %in% keepB, ]

allCorrelations <- rbindlist(list(tmpA, tmpB, tmpD))
#allCorrelations <- rbindlist(list(tmpApos, tmpAneg, tmpB, tmpDneg, tmpDpos))
##########################################
# load gene info
##########################################
geneTE<-data.table::fread("/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data/gene-TE.ID", sep = '\t', 
                          data.table = FALSE, fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
geneTE$V5<-gsub("Name=",'', geneTE$V5)
geneTE<-unlist(strsplit(geneTE$V5, split = ","))
geneTE_presence <- allCorrelations$geneID %in% geneTE
allCorrelations$geneType <- ifelse(geneTE_presence, "gene-TE", "gene")
#subset(allCorrelations, geneID=="Solyc01g111180.2") # should be gene-TE
##########################################
# get density per genePart
##########################################
allCorrelations$corrResult <- ifelse((allCorrelations$pValpearsonExpMeth>0.05 |
                                        is.na(allCorrelations$pValpearsonExpMeth)),
                                     "non-correlated", "correlated")
# Apply the function to the genePart column and store the results in a new column
allCorrelations$numeric_position <- sapply(allCorrelations$genePart, calculate_numeric_position)
# Apply the function to the genePart column
allCorrelations$gene_position <- sapply(allCorrelations$genePart, calculate_gene_position)
allCorrelations$xPos <- allCorrelations$gene_position+allCorrelations$numeric_position

values<-seq(0,max(allCorrelations$xPos), by= 1)
allCorrelations$window <- cut(allCorrelations$x, breaks = values)
allCorrelations$pos <- as.integer(gsub("\\((\\d+),(\\d+)\\]", "\\1", allCorrelations$window))

allCorrelations$ctxt<-allCorrelations$DMR
allCorrelations <- tidyr::separate(allCorrelations, ctxt, into = c("chr","dmr", "start"), sep = "_", convert = TRUE)
allCorrelations$chr<-NULL
allCorrelations$start<-NULL
allCorrelations$corrDirection <- ifelse(allCorrelations$pearsonExpMeth>0, "positive", "negative")
head(allCorrelations)
##
data.table::fwrite(allCorrelations, paste0(outDir,"06.01_correlations-per-window.bed"), 
                   quote=F, row.names=F,col.names = T,sep="\t")
gc()
####################################################################################
####################################################################################
####################################################################################
dmrOverTEs<-data.table::fread("DMRs_over_TEs.bed", sep = '\t', data.table = T, 
                              fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
dmrOverTEs
toCheckDMRs<-unique(dmrOverTEs$V4)
setDT(allCorrelations)
allCorrelations[, overTE := ifelse( DMR %in% toCheckDMRs, 'overTEs','nonTEs' ) ]

allCorrelations$corrResult <- factor(allCorrelations$corrResult)

totalDMRoverTEs <- allCorrelations %>% 
  group_by(pos, geneType,dmr,overTE) %>%
  dplyr::summarise(nDMR = length(unique(DMR))) %>%
  group_by(pos)

dmr_data <- data.table(dmr = c("C-DMR", "CG-DMR"),totalDMR = c(1794110, 293047))

totalDMRoverTEs <- left_join(totalDMRoverTEs, dmr_data, by = c("dmr"))
totalDMRoverTEs$percentage<-(totalDMRoverTEs$nDMR/totalDMRoverTEs$totalDMR)*100

ggplot(totalDMRoverTEs, aes(x = pos, y = percentage, color= overTE, fill=overTE)) +
  geom_bar(position='stack', stat='identity')+
  facet_grid( dmr ~ geneType, scales = "free") +
  scale_color_manual(values=c("white","white"))+
  scale_fill_manual(values=c("grey50","grey80"))+
  labs(x = "Position", y = "Percentage of DMRs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outPlot,"06.01_DMR_number_over_TEs-per-position.pdf"), width = 40, height = 20, units = "cm")

totalDMR <- allCorrelations %>% 
  group_by(pos, geneType,dmr) %>%
  dplyr::summarise(sumNdmr = length(unique(DMR))) %>%
  group_by(pos)

count_data <- allCorrelations %>% 
  group_by(pos, corrResult, geneType,corrDirection,dmr) %>%
  dplyr::summarise(count = length(unique(DMR))) %>%
  group_by(pos)

#### 
allCorrelations %>%
  filter(corrResult=="correlated")
resulsPerGenes <- allCorrelations %>%
  filter(pos >= 20 & pos < 24) %>%
  group_by(geneID) %>%
  dplyr::summarize(allCorr = paste(unique(corrResult), collapse = ","),
                   allDirection = paste(unique(corrDirection), collapse = ","),
                   allDMR = paste(unique(dmr), collapse = ",")) %>%
  mutate(allCorr = ifelse(allCorr %in% c("non-correlated,correlated", "correlated,non-correlated"), 'correlated', allCorr),
         allDirection = ifelse(allDirection %in% c("positive,negative", "negative,positive"), 'both', allDirection),
         allDMR = ifelse(allDMR %in% c("C-DMR,CG-DMR", "CG-DMR,C-DMR"), 'both', allDMR))

count_data <- resulsPerGenes %>% 
  group_by(allCorr) %>%
  dplyr::summarise(count = length(unique(geneID)))
#
nSigGenes <- data.frame(
  correlation = c("significative", "non-significative"),
  numberOfGenes = c(14069,18612),
  percentage= c(43.05,56.95))

ggplot(nSigGenes, aes(x = "", y = percentage, fill = correlation)) +
  geom_bar(stat = "identity", width = 1) +  
  geom_text(aes(label = paste0(round(percentage, 1), "% (", numberOfGenes, ")")), 
            position = position_stack(vjust = 0.5), color="white") +  #
  coord_polar("y", start = 0) +  # Make the chart polar
  scale_fill_manual(values = c("grey10","#46B951")) +
  theme_void()  # Remove unnecessary elements
ggsave(paste0(outPlot,"00.1_number_genes_expression-affected.pdf"), width = 40, height = 20, units = "cm")

geneCtxt <- data.frame(
  DMR = c("CG-DMR", "C-DMR","CG-/C-DMR"),
  numberOfGenes = c(1410,1083,11576),
  positive=c(226,84,243),
  negative=c(189,103,297),
  undefined=c(995,896,10936))

ggplot(geneCtxt, aes(x = "", y = numberOfGenes, fill = DMR)) +
  geom_bar(stat = "identity", width = 1) +  
  labs(title = "Number of Genes by DMR", x = "DMR", y = "Number of Genes") +
  theme_minimal() +
  geom_text(aes(label = numberOfGenes), position = position_stack(vjust = 0.5), 
            color = "black") +
  scale_fill_manual(values=c("#ffca7b","grey50","#820a86"))
ggsave(paste0(outPlot,"00.2_genes-expression-affected_ctxt.pdf"), width = 40, height = 20, units = "cm")


geneCtxt_long <- geneCtxt %>%
  tidyr::pivot_longer(cols = c("positive", "negative", "undefined"), names_to = "Direction", values_to = "Count") %>%
  group_by(DMR) %>%
  mutate(proportion = Count / sum(Count)) %>%
  ungroup()

# Create the stacked bar plot
ggplot(geneCtxt_long, aes(x = DMR, y = proportion, fill = Direction)) +
  geom_bar(stat = "identity", width = 1) +  
  labs(title = "Gene Context by DMR", x = "DMR", y = "Proportion of Genes") +
  scale_fill_manual(values = c("positive" = "#3BC000",
                               "negative" = "#F8070A", 
                               "undefined" = "#2D38D2")) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "white") +
  theme_minimal()
ggsave(paste0(outPlot,"00.3_genes-expression-affected_corrDirection.pdf"), width = 40, height = 20, units = "cm")

library(scales)
# Given values

values <- c(189, 103, 297)
values <- c(995, 886, 10937)
# Normalize the values to range from 0 to 1

# Generate a color gradient (e.g., from blue to red)
color_gradient <- colorRampPalette(c("blue", "#DF2E03"))(length(values))
values <- c(84,226, 342)
color_gradient <- colorRampPalette(c("#f7b9aa", "red"))(length(values))

normalized_values <- rescale(values)
mapped_colors <- color_gradient[rank(normalized_values)]
# Convert the mapped colors to RGB
rgb_colors <- col2rgb(mapped_colors)
rgb_codes <- apply(rgb_colors, 2, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
# Print the RGB codes
print(rgb_codes)

"#F7B9AA" "#FB5C55" "#FF0000"

toGetExpression<-allCorrelations %>%
  filter(pos>=15 & pos<24) %>%
  select(geneID,dmr,corrResult,corrDirection) %>%
  group_by(geneID) %>%
  dplyr::summarize(allCorr = paste(unique(corrResult), collapse = ","),
                   allDirection = paste(unique(corrDirection), collapse = ","),
                   allDMR = paste(unique(dmr), collapse = ",")) %>%
  mutate(allCorr = ifelse(allCorr %in% c("non-correlated,correlated", "correlated,non-correlated"), 'correlated', allCorr),
         allDirection = ifelse(allDirection %in% c("positive,negative", "negative,positive"), 'both', allDirection),
         allDMR = ifelse(allDMR %in% c("C-DMR,CG-DMR", "CG-DMR,C-DMR"), 'both', allDMR)) %>%
  filter(allCorr=="correlated") %>%
  select(geneID,allDMR,allDirection)

data.table::fwrite(toGetExpression, paste0(outDir,"genes_expression_affected_methylation.bed"), 
                   quote=F, row.names=F,col.names = T,sep="\t")
toGetExpressionNonCor<-allCorrelations %>%
  filter(pos>=15 & pos<24) %>%
  select(geneID,dmr,corrResult,corrDirection) %>%
  group_by(geneID) %>%
  dplyr::summarize(allCorr = paste(unique(corrResult), collapse = ","),
                   allDirection = paste(unique(corrDirection), collapse = ","),
                   allDMR = paste(unique(dmr), collapse = ",")) %>%
  mutate(allCorr = ifelse(allCorr %in% c("non-correlated,correlated", "correlated,non-correlated"), 'correlated', allCorr),
         allDirection = ifelse(allDirection %in% c("positive,negative", "negative,positive"), 'both', allDirection),
         allDMR = ifelse(allDMR %in% c("C-DMR,CG-DMR", "CG-DMR,C-DMR"), 'both', allDMR)) %>%
  filter(allCorr=="non-correlated") %>%
  select(geneID,allDMR,allDirection)

data.table::fwrite(toGetExpressionNonCor, paste0(outDir,"genes_expression_non-affected_methylation.bed"), 
                   quote=F, row.names=F,col.names = T,sep="\t")
###

correlatedData <- subset(count_data, corrResult=="correlated")

wide_data <- correlatedData %>%
  tidyr::pivot_wider(names_from = corrDirection, values_from = count, 
                     values_fill = list(count = 0)) %>%
  dplyr::rename(Negative = negative, Positive = positive)

wide_data <- left_join(wide_data, totalDMR, by = c("pos", "geneType", "dmr"))
wide_data <- left_join(wide_data, dmr_data, by = c("dmr"))
wide_data$percentage<-(wide_data$sumNdmr/wide_data$totalDMR)*100

wide_data <- mutate(wide_data, Negative = if_else(Negative ==0, 1, Negative),
                    Positive = if_else(Positive == 0, 1, Positive),
                    TotalAssociated = Positive + Negative)

odds_data <- wide_data %>%
  mutate(Odds_Positive = (Positive/ sumNdmr) / (TotalAssociated/sumNdmr),
         Odds_Negative = (Negative / sumNdmr) / (TotalAssociated/sumNdmr),
         OddsRatio = Odds_Positive / Odds_Negative)

odds_data$OddsRatio<-ifelse(odds_data$pos>24, 1, odds_data$OddsRatio)
ggplot(odds_data, aes(x = pos, y = OddsRatio, color = dmr, group = dmr)) +
  geom_line() +
  scale_color_manual(values=c("#820a86","#ffca7b"))+
  labs(title = "", x = "Position", y = "odd ratio",color = "Type") +
  facet_grid(.~geneType)+
  theme_bw() 
ggsave(paste0(outPlot,"01_odds_ratio.pdf"), width = 40, height = 20, units = "cm")

stacked_data <- odds_data %>%
  select(pos, geneType,dmr, Negative, Positive, percentage) %>%
  tidyr::pivot_longer(cols = c("Negative", "Positive"),
               names_to = "CountType",
               values_to = "Counts")
summary(stacked_data$Counts)

ggplot(stacked_data, aes(x = factor(pos), y = CountType, fill = Counts)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradientn(
    colours = hcl.colors(10, "OrRd", rev = TRUE),
    trans = "log10",
    limits = c(1, max(stacked_data$Counts, na.rm = TRUE)),
    breaks = 10^seq(0, ceiling(log10(max(stacked_data$Counts))), by = 1)
  ) +
  labs(
    title = "Number of Associated DMRs per Position",
    fill = "Log of Counts"
  ) +
  geom_text(aes(label = round(Counts,1)), color = "white", size = 2, angle = 90) +  # Add this line to include text labels
  coord_fixed(ratio = 4) +
  facet_grid(dmr~geneType)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5))
ggsave(paste0(outPlot,"02_heatmap_associated_DMR-per-positions.pdf"), width = 40, height = 20, units = "cm")


#
correlatedData <- subset(allCorrelations, corrResult=="correlated")
meanCorr <- correlatedData %>% 
  group_by(pos,geneType,dmr,overTE,corrDirection) %>%
  summarise(meanCorrelation = mean(pearsonExpMeth)) %>%
  group_by(pos)
meanCorr$label<-paste0(meanCorr$overTE,':',meanCorr$corrDirection)

meanCorr$meanCorrelation<-ifelse(meanCorr$pos>24, 0, meanCorr$meanCorrelation)
ggplot(meanCorr, aes(x = factor(pos), y = label, fill = meanCorrelation)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradientn(
    colours = hcl.colors(10, "RdYlGn"),
    limits=c(-0.6,0.6)) +
  labs(title = "Number of Associated DMRs per Position",
    fill = "Log of Counts") +
  geom_text(aes(label = round(meanCorrelation,1)), color = "white", size = 2, angle = 90) +  # Add this line to include text labels
  coord_fixed(ratio = 4) +
  facet_grid(dmr~geneType)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5))
ggsave(paste0(outPlot,"03.a_heatmap_mean-correlation_associated_DMR-per-positions.pdf"), width = 40, height = 20, units = "cm")

meanCorr <- correlatedData %>% 
  group_by(pos,geneType,dmr,corrDirection) %>%
  summarise(meanCorrelation = mean(pearsonExpMeth)) %>%
  group_by(pos)

meanCorr$meanCorrelation<-ifelse(meanCorr$pos>24, 0, meanCorr$meanCorrelation)

ggplot(meanCorr, aes(x = factor(pos), y = corrDirection, fill = meanCorrelation)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_fill_gradientn(
    colours = hcl.colors(10, "RdYlGn"),
    limits=c(-0.6,0.6)) +
  labs(title = "Number of Associated DMRs per Position",
       fill = "Log of Counts") +
  geom_text(aes(label = round(meanCorrelation,1)), color = "white", size = 2, angle = 90) +  # Add this line to include text labels
  coord_fixed(ratio = 4) +
  facet_grid(dmr~geneType)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5))
ggsave(paste0(outPlot,"03.b_heatmap_mean-correlation_associated_DMR-per-positions_woTEs.pdf"), width = 40, height = 20, units = "cm")
