suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
})
# Input and output directories
inDir<-"/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.0_get-DMR/ab_predigree/"
outDir<-"/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree/"

input<-list.files(path=inDir, pattern = '.bed', full.names = FALSE)
df<-c()
#pairWiseStat<-c()
for (i in input){
  DMR<-unlist(strsplit(i, '.', fixed = T))[2]
  toGen<-unlist(strsplit(i, '.', fixed = T))[1]
  
  chr<-unlist(strsplit(toGen, '_', fixed = T))[1]
  
  g1<-unlist(strsplit(toGen, '_', fixed = T))[2]
  g1<-unlist(strsplit(g1, '-', fixed = T))[2]
  g2<-unlist(strsplit(toGen, '_', fixed = T))[3]
  g2<-unlist(strsplit(g2, '-', fixed = T))[2]
  
  s1<-unlist(strsplit(toGen, '_', fixed = T))[2]
  s2<-unlist(strsplit(toGen, '_', fixed = T))[3]
  s1<-unlist(strsplit(s1, '-', fixed = T))[1]
  s2<-unlist(strsplit(s2, '-', fixed = T))[1]
  
  if(g1==g2){
    comparison<-paste0("intra-", g1)

  }else {
    comparison<-"inter-Gen"
  }
  tmp <- tryCatch({
    fread(paste0(inDir, i), sep = '\t', header = F, 
          data.table = T, fill = TRUE, na.strings = c(""), nThread = 20)
  }, warning = function(w) {
    message(sprintf("Warning while reading file %s: %s", i, conditionMessage(w)))
    NULL  # Return NULL if there's a warning
  }, error = function(e) {
    message(sprintf("Error while reading file %s: %s", i, conditionMessage(e)))
    NULL  # Return NULL if there's an error
  })
  if (is.null(tmp) ) {
    next
  }
  #Get the biary table
  binTmp<-tmp[,c(1,2,3,12,13,14)]
  binTmp$sample1<-tmp$V6/tmp$V5
  binTmp$sample2<-tmp$V9/tmp$V8
  nWin<-nrow(binTmp) # number of MethRegions
  binTmp<- binTmp %>%
    filter(V13<=0.05) %>%
    #rename_with(~ c(s1, s2), c("sample1", "sample2")) %>%
    rowwise() %>%
    mutate(
      median_meth = median(c_across(c("sample1", "sample2")), na.rm = TRUE), 
      across(all_of(c("sample1", "sample2")), ~ ifelse(. >= median_meth, 1, 0), .names = "bin_{.col}")
    ) %>%
    select(-median_meth) %>%
    ungroup()
  # if(dim(binTmp)[1]==0){
  #   nDMR<-NA
  #   transition_counts <- tibble(
  #     from_0_to_1 = NA,
  #     from_1_to_0 = NA
  #   )
  # }else{
  #   nDMR<-nrow(binTmp)
  #   # Count the number of transitions
  #   # binTmp <- binTmp %>%
  #   #   summarise(
  #   #     from_0_to_1 = sum(`bin_sample1` == 0 & `bin_sample2` == 1, na.rm = TRUE),
  #   #     from_1_to_0 = sum(`bin_sample1` == 1 & `bin_sample2` == 0, na.rm = TRUE)
  #   #   )
  # }
  samplePair<-paste0(s1,':',g1,'-',s2,':',g2)
  #tmpStat<-cbind.data.frame(DMR, chr,comparison, samplePair,nWin,nDMR,transition_counts)
  #pairWiseStat<-rbind.data.frame(pairWiseStat, tmpStat)
  # Skip to next iteration if the file is NULL
  if (nrow(binTmp) == 0) {
      next
  }
  binTmp$samplePair<-samplePair
  #bin_name <- paste0(outDir, DMR, "_",chr,'_', comparison, "_", samplePair, ".forBinary")
  #data.table::fwrite(binTmp, file = bin_name, quote = FALSE, 
  #                   row.names = FALSE, col.names = T, sep = "\t")
  binTmp$DMR<-DMR
  binTmp$comparison<-comparison
  binTmp$sample<-paste0(s1,'-',s2)
  df<-rbind.data.frame(df,binTmp)
}
##### get the loci
df$pos<-paste0(gsub("SL2.50",'',df$V1),"-",df$V2,"-",df$V3)
  setorder(df, pos)

########################################################################
###### get coordinates per comparison
########################################################################
{
  summDT <- df %>%
    group_by(pos) %>%
    summarise(DMRs = paste(unique(DMR), collapse = ";"),  
              nPairSamples = paste(unique(sample), collapse = ";"), 
              meanMethDiff = mean(V14, na.rm = TRUE),  # Calculate mean of V14
               .groups = "drop") %>%
    separate(pos, into = c("chr", "start", "end"), sep = "-", convert = TRUE)%>%
    mutate(sample_count = lengths(strsplit(nPairSamples, ";")),
           DMRs = ifelse(DMRs == "C-DMR;CG-DMR","C-DMR",
                         ifelse(DMRs == "CG-DMR;C-DMR",
                                "C-DMR", DMRs)))
  for(d in unique(summDT$DMRs)){
    tmp<-summDT%>%
      filter(DMRs==d)
    tmp$chr<-paste0("SL2.50", tmp$chr)
    data.table::fwrite(tmp, file = paste0(outDir, d, "_all.positions"), quote = FALSE, 
                       row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  
# Summarize and group by pos and comparison
summDT <- df %>%
  group_by(pos, comparison) %>%
  summarise( DMRs = paste(unique(DMR), collapse = ";"),  
    nPairSamples = paste(unique(sample), collapse = ";"), 
    meanMethDiff = mean(V14, na.rm = TRUE),  # Calculate mean of V14
    .groups = "drop") %>%
  separate(pos, into = c("chr", "start", "end"), sep = "-", convert = TRUE)%>%
  mutate(sample_count = lengths(strsplit(nPairSamples, ";")),
         DMRs = ifelse(DMRs == "C-DMR;CG-DMR","C-DMR",
                       ifelse(DMRs == "CG-DMR;C-DMR",
                              "C-DMR", DMRs)))
summDT$chr<-paste0("SL2.50", summDT$chr)
# Define the combinations of DMR types and comparisons
# save loci for each comparison
dmr_types <- unique(summDT$DMRs)
comparisons <- unique(summDT$comparison)
for (dmr in dmr_types) {
  for (comp in comparisons) {
    tmp <- summDT %>%
      filter(DMRs == dmr, comparison == comp) %>%
      #select(chr, start, end) %>%
      arrange(chr, start, end)
    file_name <- paste0(outDir, dmr, "_", gsub("-", "", comp), ".positions")
    data.table::fwrite(tmp, file = file_name, quote = FALSE, 
                       row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}
  
# Count the number of positions for each comparison and sample_count
pedResult<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/af_pedigree/"
pedPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/plots/"
result <- summDT %>%
  group_by(comparison, sample_count, DMRs) %>%
  summarise(num_positions = n(), .groups = "drop")
data.table::fwrite(result, file = paste0(pedResult,"05.1_number_variable_positions.tsv"), quote = FALSE, 
                   row.names = FALSE, col.names = T, sep = "\t")
data.table::fwrite(df, file = paste0(pedResult,"05.0_DMRs-pedigree.tsv"), quote = FALSE, 
                   row.names = FALSE, col.names = T, sep = "\t")

}
#############################
#
#
pairWiseStatSumm<-df %>%
    group_by(DMR, comparison,samplePair) %>%
      summarise(
        nDMR = n(),
        from_0_to_1 = sum(`bin_sample1` == 0 & `bin_sample2` == 1, na.rm = TRUE),
        from_1_to_0 = sum(`bin_sample1` == 1 & `bin_sample2` == 0, na.rm = TRUE),
        .groups = 'drop')

pairWiseStatSumm$comparison <- factor(pairWiseStatSumm$comparison, levels = c("intra-G0", "intra-G5", "inter-Gen"))
plot <- ggplot(pairWiseStatSumm, aes(x = comparison, y = nDMR, fill = comparison, color = comparison)) +
  geom_boxplot(aes(color = comparison),outlier.shape = NA, alpha = 0.6) +  # Boxplot without default outliers
  geom_jitter(aes(color = comparison), width = 0.2, size = 2, alpha = 0.6) +  # Add dots with jitter
  facet_grid(DMR ~ ., scales = "free_y") +
  scale_fill_manual(values = c("#aeaeae", "#727172", "#232323")) +
  scale_color_manual(values = c("#aeaeae", "#727172", "#232323")) +
  labs(
    title = "# of DMR by pairwise comparison",
    x = "Comparison",
    y = "# of DMR"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(paste0(pedPlot, "05.1_DMR-number_between-generations.pdf"), plot = plot, width = 8, height = 6)



custom_order <- c("C3:G0", "C6:G0", "C9:G0", "C1:G5", "C2:G5", "C3:G5")
contingency_table <- pairWiseStatSumm %>%
  mutate(
    group1_raw = sub("-.*", "", samplePair),
    group2_raw = sub(".*-", "", samplePair),
    group1 = factor(
      ifelse(match(group1_raw, custom_order) <= match(group2_raw, custom_order), group1_raw, group2_raw),
      levels = custom_order
    ),
    group2 = factor(
      ifelse(match(group1_raw, custom_order) > match(group2_raw, custom_order), group1_raw, group2_raw),
      levels = custom_order
    )
  ) %>%
  group_by(DMR, group1, group2, comparison) %>%
  summarise(total_nDMR = sum(nDMR, na.rm = TRUE), .groups = "drop")


print(contingency_table, n = 30)
heatmap_long <- contingency_table %>%
  mutate(comparison = factor(comparison, levels = c("intra-G0", "inter-Gen", "intra-G5"))) %>%
  filter(DMR=="C-DMR")
plot <- ggplot(heatmap_long, aes(x = group1, y = group2, fill = total_nDMR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = total_nDMR), color = "black", size = 3) + # Add numbers in the boxes
  scale_fill_gradient(low = "yellow", high = "red", name = "# of nDMR") +
  labs(title = "Heatmap of DMR Counts by Sample Comparison", x = "Group 1", y = "Group 2") +
  theme_minimal() +
  facet_wrap(~ DMR, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(pedPlot, "05.2_C-DMR_nDMR_generation.pdf"), plot = plot, width = 8, height = 6)

heatmap_long <- contingency_table %>%
  mutate(comparison = factor(comparison, levels = c("intra-G0", "inter-Gen", "intra-G5"))) %>%
  filter(DMR=="CG-DMR")
plot <- ggplot(heatmap_long, aes(x = group1, y = group2, fill = total_nDMR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = total_nDMR), color = "black", size = 3) + # Add numbers in the boxes
  scale_fill_gradient(low = "yellow", high = "red", name = "# of nDMR") +
  labs(title = "Heatmap of DMR Counts by Sample Comparison", x = "Group 1", y = "Group 2") +
  theme_minimal() +
  facet_wrap(~ DMR, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(pedPlot, "05.2_CG-DMR_nDMR_generation.pdf"), plot = plot, width = 8, height = 6)


##########################################################################################
# get number of changes
##########################################################################################
methBin<-df %>%
  select(pos,DMR,comparison, samplePair, sample1,sample2, bin_sample1, bin_sample2)
methBin <- methBin %>%
  filter(comparison=="inter-Gen") %>%
  mutate(
    # Split samplePair into two parts
    sample1_gen = sub(".*:(G[0-9]+).*", "\\1", samplePair),
    sample2_gen = sub(".*:(G[0-9]+)$", "\\1", samplePair),
    
    # Determine ancestral/descendant state
    state_sample1 = ifelse(sample1_gen == "G0", "ancestral", "descendant"),
    state_sample2 = ifelse(sample2_gen == "G0", "ancestral", "descendant"),
    
    # Combine with bin_sample columns to indicate state
    ancestral = ifelse(state_sample1 == "ancestral", bin_sample1, bin_sample2),
    descendant = ifelse(state_sample1 == "descendant", bin_sample1, bin_sample2)
  )

# Step 1: Count the number of discordant samples per position
discordant_counts <- methBin %>%
  group_by(pos, DMR) %>%
  summarise(
    num_discordant = sum(ancestral != descendant),
    num_meth = sum(descendant>0),
    num_unmeth = sum(descendant<1))

# Step 2: Count the number of positions per number of discordant samples
discordant_summary <- discordant_counts %>%
  group_by(DMR,num_discordant,num_meth, num_unmeth) %>%
  summarise(
    num_positions = n(), # Count positions for each discordant count
    .groups = "drop"
  )


discordant_long <- discordant_counts %>%
  pivot_longer(
    cols = c(num_meth, num_unmeth),
    names_to = "state",          # Create a new column for state (methylated/unmethylated)
    values_to = "count"          # Value column for counts
  )

# Step 2: Summarize counts by num_discordant and state
discordant_summary <- discordant_long %>%
  group_by(DMR,num_discordant, state) %>%
  summarise(
    total_positions = sum(count), # Aggregate counts for each state
    .groups = "drop"
  )

# Create the stacked bar plot
ggplot(discordant_summary, aes(x = as.factor(num_discordant), y = total_positions, fill = state)) +
  geom_bar(stat = "identity", position = "stack", color = "black") + # Stacked bar chart
  labs(
    x = "# of G5 plants with discordant meth state with G0 plants",
    y = "Number of Positions",
    fill = "State",
    title = "Methylated and Unmethylated DMRs"
  ) +
  facet_grid(DMR~., scales = "free")+
  scale_fill_manual(
    values = c("num_meth" = "black", "num_unmeth" = "white"),
    labels = c("Methylated", "Unmethylated")
  ) +
  theme_minimal()
ggsave(paste0(pedPlot, "05.3_binary-state.pdf"), width = 8, height = 6)




