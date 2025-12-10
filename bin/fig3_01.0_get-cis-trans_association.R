setwd("/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/ba_snp-blocks/cb_result")

outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/result"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/plots"
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(eulerr)
  library(dplyr)
  library(ggplot2)
})
# Read data
input<-list.files(pattern = ".mQTL", full.names = FALSE)
dt <- lapply(input, function(file_path) {
  out <- fread(file_path, sep = '\t', header = F, fill = TRUE, na.strings = "",nThread = 40)
  out[, diff_V2 := c(0, diff(V2)), by = V1] # get consecutive distance of snps for each chr
  out[, groups := cumsum(diff_V2 > 500000), by = V1] # get the group if they are 50000 bp away
  # Merge rows within each group
  result <- out[, .(min_V2 = min(V2), # Keep the minimum value of SNP start
                   max_V3 = max(V3), # Keep the maximum value of SNP end
                   mean_V4 = mean(V4), # Average of V4
                   mean_V5 = mean(V5), # Average of V5
                   min_V6 = min(V6), # Keep the minimum value of V6 PVALUE
                   nSNPs= .N,
                   V7 = toString(unique(V7))), # Collapse unique values of V7 by comma
               by = .(V1, groups)]  #use both
  # Remove the 'groups' column
  #result[, groups := NULL]
  # Remove the 'diff_V2' column
  #result[, diff_V2 := NULL]
  result <- result[, .(V1, min_V2, max_V3, mean_V4, mean_V5, min_V6,nSNPs, V7 = unlist(strsplit(V7, ","))), by = 1:nrow(result)]
  result[, c("chrDMR", "DMR", "posDMR") := tstrsplit(V7, "_", fixed = TRUE)]
  result$chrDMR<-as.numeric(gsub("ch",'',result$chrDMR))
  result$posDMR<-as.numeric(result$posDMR)
  result[, lociType := ifelse(chrDMR == V1 & posDMR <= max_V3 & posDMR >= min_V2, "cis", "trans")]
  result[, lociPos := ifelse(chrDMR == V1, "same-chr", "diff-chr")]
  return(result)
})
dt <- rbindlist(dt)
dt[, nrow := NULL]
colnames(dt)<-c('chrSNP','start-loci','end_loci','meanBeta','meanSDbeta','minPvalue',
                'nSNPs','DMRpos','chrDMR','DMR','startDMR','lociType','lociPosition')

data.table::fwrite(dt, file= paste0(outDir,'01.0_cis-trans_sig-associations.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")
# do the lollypops
# Adjust the categorization to make sure all relevant categories are covered
total_counts <- dt %>%
  group_by(DMR) %>%
  summarise(total_DMRs = n_distinct(DMRpos), .groups = 'drop')

mQTL_number<- dt%>%
  group_by(DMRpos, DMR) %>%
  dplyr::summarize(mQTL_number = n()) %>%
  mutate(class = case_when(
    mQTL_number == 1 ~ "1",
    mQTL_number > 1 & mQTL_number <= 10 ~ "2-10",
    mQTL_number > 10 & mQTL_number <= 20 ~ "11-20",
    mQTL_number > 20 & mQTL_number <= 30 ~ "21-30",
    mQTL_number > 30 & mQTL_number <= 40 ~ "31-40",
    mQTL_number > 40 & mQTL_number <= 50 ~ "41-50",
    mQTL_number > 50 & mQTL_number <= 100 ~ "51-100",
    mQTL_number > 100 ~ ">100"
  ))

class_counts <- mQTL_number %>%
  group_by(DMR, class) %>%
  summarise(count = n(), .groups = 'drop')

proportion_counts <- class_counts %>%
  left_join(total_counts, by = c("DMR")) %>%
  mutate(proportion = count / total_DMRs)

#
ggplot(proportion_counts, aes(x = proportion, y = class, fill = DMR,color = DMR)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1.0), alpha = 0.7, width = 0.9) +
  labs(title = "Proportion of DMRs with mQTL in Each Category",
       x = "Proportion of DMRs",
       y = "Numbers of mQTL") +
  scale_fill_manual(values = c("#820a86", "#ffca7b")) +
  scale_color_manual(values = c("#820a86", "#ffca7b")) +
  theme_minimal() +
  scale_y_discrete(limits = c("1", "2-10", "11-20", "21-30", "31-40", "41-50", "51-100", ">100")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "top")
ggsave(paste0(outPlot,"01.0_proportion_DMRs-per-mQTLnumber.pdf"), width = 12, height = 8)

# # Create the lollipop plot
# dodge_width <- 0.5
# # Create the lollipop plot with dodged positions to avoid stacking
# ggplot(proportion_counts, aes(x = proportion, y = class, color = DMR)) +
#   geom_segment(aes(yend = class, x = 0, xend = count), size = 1, position = position_dodge(width = dodge_width)) +
#   geom_point(size = 3, position = position_dodge(width = dodge_width)) +
#   scale_y_discrete(limits = c("1", "2-10", "11-20", "21-30", "31-40", "41-50", "51-100", ">100")) +  # Ensure y-axis categories are ordered
#   labs(title = "Lollipop Plot of mQTL Number Distribution per DMR",
#        x = "Proportion of DMRs",
#        y = "Numbers of mQTL") +
#   scale_color_manual(values = c("#820a86", "#ffca7b")) +
#   theme_minimal() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1))
# ggsave(paste0(outDir,"lollypop_counts.pdf"),  width = 40, height = 20, units = "cm")

summ <- dt %>%
  group_by(DMRpos, DMR) %>%
  mutate(DMRpos = gsub(" ", "", DMRpos)) %>%
  dplyr::summarise(type = paste0(unique(lociType),collapse = ","),
                   nSNPs = n(),
                   chrs = paste0(unique(chrSNP), collapse = ",")) %>%
  mutate(type = if_else(grepl("cis", type) & grepl("trans", type), "cis,trans", type))

data.table::fwrite(summ, file= paste0(outDir,'01.1_cis-sig-associations_summary.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")
## make a venn
summ %>%
  group_by(DMR, type) %>%
  dplyr::summarize(count = n())
# A tibble: 6 × 3
# Groups:   DMR [2]
# DMR    type      count
# <chr>  <chr>     <int>
# 1 C-DMR  cis          80
# 2 C-DMR  cis,trans   455
# 3 C-DMR  trans     10460
# 4 CG-DMR cis          43
# 5 CG-DMR cis,trans   199
# 6 CG-DMR trans      4441
# getVennCGdmr <- euler(c("cis" = 43, "trans" = 4441, "cis&trans" = 199),shape = "ellipse")
# pdf(paste0(outDir, "01.1_CG-DMR_cis-trans_venn.pdf"))
# plot(getVennCGdmr,quantities = list(type = c("percent", "counts"), font = 3))
# dev.off()

pieData<- data.frame(
  DMR = rep(c('C-DMR', 'CG-DMR'), each=2),
  category = c("cis", "trans"),
  count = c(535, 10460,
            242, 4441))
pieData$category <- factor(pieData$category,
                           levels = c("cis", "trans"))

# Calculate total and proportion
pieData <- pieData %>%
  group_by(DMR) %>%
  mutate(total = sum(count),
         proportion = count / total * 100)

ggplot(pieData, aes(x = "", y = proportion, fill = category)) +
  geom_bar(width = 1, stat = "identity", alpha=0.7) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white") +
  labs(title = "Proportion of DMR Categories", fill = "Category") +
  #scale_fill_manual(values = c("#3bb6c4","grey70","grey50","darkred",'#ffa56a'))+
  scale_fill_manual(values = c("grey20","grey50"))+
  theme_void() +
  facet_grid(~DMR) +
  theme(legend.position = "right")
ggsave(paste0(outPlot, "01.1_DMRs-association_cis-trans_piePlot.pdf"), width = 12, height = 8)
#
# A tibble: 6 × 4
# DMRtype associationType  count percentage
# <fct>   <fct>            <int>      <dbl>
# 1 C-DMR   NonSig          239472       56.1
# 2 C-DMR   NonSig_by_GIF   176211       41.3
# 3 C-DMR   sig_wo_GIF       10995        2.6
# 4 CG-DMR  NonSig           53628       49.8
# 5 CG-DMR  NonSig_by_GIF    49460       45.9
# 6 CG-DMR  sig_wo_GIF        4683        4.3

pieData<- data.frame(
  DMR = rep(c('C-DMR', 'CG-DMR'), each=4),
  category = c("unlinked", "linked with GI","cis", "trans"),
  count = c(239664, 176211, 535, 10460,
            53531,49560, 242, 4441))
pieData$category <- factor(pieData$category, levels = c("unlinked", "linked with GI", "cis", "trans"))

# Calculate total and proportion
pieData <- pieData %>%
  group_by(DMR) %>%
  mutate(total = sum(count),
         proportion = count / total * 100)
CHANGE COLORS
ggplot(pieData, aes(x = "", y = proportion, fill = category)) +
  geom_bar(width = 1, stat = "identity", alpha=0.7) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white") +
  labs(title = "Proportion of DMR Categories", fill = "Category") +
  scale_fill_manual(values = c("#3bb6c4","grey70","grey50","darkred",'#ffa56a'))+
  theme_void() +
  facet_grid(~DMR) +
  theme(legend.position = "right")
ggsave(paste0(outPlot, "01.2_DMRs-association_summary_piePlot.pdf"), width = 12, height = 8)




