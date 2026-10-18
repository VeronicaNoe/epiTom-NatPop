suppressPackageStartupMessages({
  library(ggridges)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gplots)
  library(RColorBrewer)
  library(rlang)
  library(pheatmap)
  library(tibble)
  
})

# Define Directories
inDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/bc_predicted-values/different-clustering/"
outDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/results/"
outPlot <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/plots/"

########################## ##########################
# A Load and Process Data
####################################################
{
load_data <- function(file) {
  fread(file, sep = '\t', data.table = T, header = T, fill = TRUE, na.string = "NA", nThread = 20)
}

# Load Master Table
masterTable <- load_data(paste0(outDir, '08.1_master_table.tsv')) %>%
  filter(kinship == "DMR-SNP")

traitGWAS <- load_data(paste0(outDir, '02.3_trait_H2.tsv')) %>%
  filter(kin == "DMR-SNP") %>%
  select(meta, h2, pseudo) %>%
  distinct() %>%
  arrange(meta)

metaInfo <- load_data(paste0(outDir, '02.3_trait_H2.tsv')) %>%
  filter(kin == "DMR-SNP") %>%
  select(-h2, -pseudo, -kin) %>%
  arrange(meta) %>%
  group_by(meta) %>%
  summarise(association = ifelse("sig" %in% association, "sig",
                                 ifelse("nonSig-GIF" %in% association, "nonSig-GIF", "NonSig"))) %>%
  left_join(traitGWAS, by = c("meta" = "meta")) %>%
  arrange(meta)

metaName <- load_data("/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype/raw/metabolites_name.csv")
metaInfo <- metaName %>% left_join(metaInfo, by = c("Name" = "meta"))

########################################
# Load Prediction Data
########################################
#dataSet<-c("PCA-balance", "PCA-number-balance","random-number-balance","unbalance")
#PCA-number-balance: 80k-both  
#PCA-balance: 120k-snp-80k-dmr.tsv  
#random-number-balance: wo-na-1313-loci.tsv   
#unbalance: snp-wo-na-all-dmrs.tsv  #original
dataSet<-c("unbalance")
for( d in dataSet){
  input <- list.files(path = paste0(inDir,d), pattern = "LOCO", full.names = FALSE)
  traitP <- do.call(rbind, lapply(input, function(i) load_data(paste0(inDir,d,"/", i))))
  datasetType<-unlist(strsplit(input[1],'.', fixed = T))[1]
  datasetType<-unlist(strsplit(datasetType,'_', fixed = T))[4]
  traitP$Cluster_Type<-NULL
  traitP$Group<-d
  traitP$Dataset_Type<-datasetType
  # write.table(traitP, file = paste0(outDir, "11.0_genomic-prediction_model_",d,".tsv"), 
  #             sep = "\t", quote = FALSE, row.names = FALSE)
  # # Compute Mean Metrics
  merged_data <- traitP %>%
    group_by(Trait, Model, Marker, Dataset_Type) %>%
    summarise(avgPA = mean(PA, na.rm = TRUE),
              avgMSPE = mean(MSPE, na.rm = TRUE),
              avgMAPE = mean(MAPE, na.rm = TRUE)) %>%
    mutate(avgPA = ifelse(avgPA < 0, 0, avgPA)) %>%
    arrange(Trait, Model, Marker, Dataset_Type) %>%
    left_join(metaInfo, by = c("Trait" = "Name"))
  
  # Categorize H2
  merged_data <- merged_data %>%
    mutate(H2_group = case_when(
      h2 <= 0.25 ~ "<= 0.25",
      h2 > 0.25 & h2 <= 0.75 ~ "> 0.25 & <= 0.75",
      h2 > 0.75 ~ "> 0.75"
    ))
  
  # Compute Mean Across Models
  mean_avg <- merged_data %>%
    group_by(Marker, Model) %>%
    summarize(mean_avgPA = mean(avgPA, na.rm = TRUE),
              mean_avgMSPE = mean(avgMSPE, na.rm = TRUE),
              mean_avgMAPE = mean(avgMAPE, na.rm = TRUE))
  
  # Compute T-test Results
  test_results <- merged_data %>%
    group_by(Model) %>%
    summarise(p_value_ttest = t.test(avgPA ~ Marker)$p.value)
  
  ########################## ##########################
  # Generalized Function for Kruskal-Wallis Test and Boxplots
  ####################################################
  run_kruskal_plot <- function(data, metric, title, y_label, file_name) {
    metric_sym <- sym(metric)  # Convert string to symbol
    
    # Run Kruskal-Wallis test
    ktest_results <- data %>%
      group_by(Marker) %>%
      summarise(p_value = kruskal.test(!!metric_sym ~ Model)$p.value, .groups = "drop")
    
    # Create a p-value dataframe for annotation
    pvalue_df <- data %>%
      group_by(Marker) %>%
      summarise(
        x = 1,  
        y = max(log10(abs(!!metric_sym)), na.rm = TRUE) * 1.05,  # Use symbol with !!
        .groups = "drop"
      ) %>%
      left_join(ktest_results, by = "Marker") %>%
      mutate(label = paste0("p-value: ", formatC(p_value, format = "e", digits = 2)))  # Fix p-value label
    
    # Create boxplot with Kruskal-Wallis test annotations
    plot <- ggplot(data, aes(x = Model, y = log10(abs(!!metric_sym)), color = Marker, fill = Marker)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("#c06a80", "#999999")) +
      scale_fill_manual(values = c("#c06a80", "#999999")) +
      facet_grid(Marker ~ ., scales = 'free') +
      labs(title = title, x = "Model", y = y_label) +
      theme_minimal() +
      geom_text(data = pvalue_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0.5, size = 5)
    
    # Save the plot
    ggsave(paste0(outPlot, file_name), plot = plot, width = 40, height = 20, units = "cm")
  }
  # Run Kruskal-Wallis tests and plot results for MSPE and MAPE
  run_kruskal_plot(merged_data, "avgMSPE", "Boxplot of MSPE by Marker and Model", "log10 Mean Squared Prediction Error (MSPE)", paste0("06.1_MSPE-boxplot_",d,".pdf"))
  run_kruskal_plot(merged_data, "avgMAPE", "Boxplot of MAPE by Marker and Model", "Mean Squared Prediction Error (MAPE)", paste0("06.1_MAPE-boxplot_",d,".pdf"))
  
  # Density Plot for Predictive Ability (PA)
  ggplot(merged_data, aes(x = avgPA, color = Marker, fill = Marker)) +
      geom_density(alpha = 0.1) +
      scale_color_manual(values = c("#c06a80", "#999999")) +
      scale_fill_manual(values = c("#c06a80", "#999999")) +
      xlim(0,1) +
      facet_grid(Model ~ Marker, scales = 'free') +
      labs(title = "Density Plot of PA by Marker and Model",
           x = "Predictive Ability (PA)", 
           y = "Density") +
      geom_vline(data = mean_avg, aes(xintercept = mean_avgPA, color = Marker), 
                 linetype = "dashed", size = 0.8) +
      theme_minimal() +
      geom_text(data = test_results, aes(x = 0.5, y = Inf, 
                                         label = paste0("p-value: ", formatC(p_value_ttest, 
                                                                             format = "e", digits = 2))),
                hjust = 1.1, vjust = 2, size = 4, inherit.aes = FALSE)
    ggsave(paste0(outPlot, "06.0_PA-density_",d,".pdf"), width = 40, height = 40, units = "cm")
    rr_reml_data <- merged_data %>%
      filter(Model == "RR-REML") %>%
      select(Trait, Marker, avgPA) %>%
      pivot_wider(names_from = Marker, values_from = avgPA) %>%
      rename(SNP_PA = SNP, DMR_PA = DMR)
    # Scatter plot
    # Compute mean values for each Marker
    mean_values <- merged_data %>%
      filter(Model == "RR-REML") %>%
      group_by(Marker) %>%
      summarise(SNP_Mean = mean(avgPA[Marker == "SNP"], na.rm = TRUE),
                DMR_Mean = mean(avgPA[Marker == "DMR"], na.rm = TRUE))
    
    # Scatter plot with vertical & horizontal mean lines
    ggplot(rr_reml_data, aes(x = SNP_PA, y = DMR_PA)) +
      geom_point(color = "black", size = 3, alpha = 0.7) +  
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  
      geom_vline(data = mean_values, aes(xintercept = SNP_Mean), 
                 linetype = "dashed", size = 0.8, color = "#999999") +  
      geom_hline(data = mean_values, aes(yintercept = DMR_Mean), 
                 linetype = "dashed", size = 0.8, color = "#c06a80") +  
      xlim(0, 1) +  # Set x-axis limits
      ylim(0, 1) +  # Set y-axis limits
      theme_minimal() +
      labs(title = "Predictive Ability (PA) Comparison: SNP vs. DMR (RR-REML)",
           x = "SNP-based PA",
           y = "DMR-based PA") +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0(outPlot, "06.0_PA_SNP_vs_DMR_RR_REML_means_lines.pdf"), width = 40, height = 40, units = "cm")
    
}
}
###########################################################
#### B get SNP-PA DMR-PA differentce per trait
##############################################################################
input <- list.files(path = paste0(outDir), pattern = "11.0_genomic-prediction_model_", full.names = FALSE)
allPA<-c()
for( i in input){
  getPA<-fread(paste0(outDir,i), sep = '\t', data.table = T, header = T, fill = TRUE, na.string = "NA", nThread = 20)
  allPA<-rbind(allPA, getPA)
}

# Compute difference (SNP PA - DMR PA) for each trait across datasets
pa_difference <- allPA %>%
  group_by(Trait, Group, Model) %>%
  summarise(avgSNP=mean(PA[Marker == "SNP"], na.rm=T),
            avgDMR=mean(PA[Marker == "DMR"], na.rm=T),
    Diff_PA = mean(PA[Marker == "SNP"], na.rm=T) - mean(PA[Marker == "DMR"], na.rm=T), .groups = "drop")


ggplot(pa_difference, aes(x = Diff_PA, color = Group, fill = Group)) +
  geom_density(alpha = 0.1) +
  #scale_color_manual(values = c("#c06a80", "#999999")) +
  #scale_fill_manual(values = c("#c06a80", "#999999")) +
  facet_grid(Model ~ ., scales = 'free') +
  labs(title = "Density Plot of PA by Marker and Model",
       x = "diff SNPs-PA DMRs-PA", 
       y = "Density") +
  theme_minimal() 
ggsave(paste0(outPlot, "06.2_diff-PA-density.pdf"), width = 40, height = 20, units = "cm")
# 
# Define color palette
color_palette <- colorRampPalette(c("#313695", "white", "#a50026"))(100)
# Loop through each model to generate heatmaps separately
for (mod in unique(pa_difference$Model)) {
  # Subset data for the current model
  model_data <- pa_difference %>%
    filter(Model == mod) %>%
    select(Trait, Group, Diff_PA) %>%
    pivot_wider(names_from = Group, values_from = Diff_PA) %>%
    column_to_rownames(var = "Trait") %>%
    as.matrix()
  
  # Handle missing values (Option: Replace NA with column median)
  for (col in 1:ncol(model_data)) {
    model_data[is.na(model_data[, col]), col] <- 0
  }
  
  # Generate heatmap for the current model
  pdf(file = paste0(outPlot, "06.3_PA_Difference_Heatmap_", mod, ".pdf"), width = 10, height = 8)
  heatmap.2(
    model_data,
    col = color_palette,
    trace = "none",
    scale = "none",
    dendrogram = "row",
    margins = c(8, 10),
    key = TRUE,
    key.title = "PA Difference",
    key.xlab = "SNP - DMR",
    density.info = "none",
    main = paste("PA Differences for", mod)
  )
  dev.off()
}
# Save the plot

formatted_table <- pa_difference %>%
  pivot_wider(
    names_from = Group,
    values_from = c(avgSNP, avgDMR, Diff_PA),
    names_glue = "{Group}_{.value}"
  ) %>%
  arrange(Trait, Model)

# Rename columns for better readability
formatted_table <- formatted_table %>%
  rename(
    Trait = Trait,
    Model = Model
  )
write.table(formatted_table, file = paste0(outDir, "12.0_dataset-comparison.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

######################################################################
##### C Characterize Traits per class
######################################################################
filtered_data <- merged_data %>%
  filter(Model == "RR-REML") %>%
  mutate(Class = factor(Class, levels = c(setdiff(unique(Class), "Others"), "Others")))

class_splits <- split(filtered_data, filtered_data$Class)
for (class_name in names(class_splits)) {
  class_data <- class_splits[[class_name]]
  # Process the data for heatmap generation
  heatmap_data <- class_data %>%
    mutate(Model_Marker = paste0(Model, "_", Marker)) %>%
    ungroup() %>%
    select(Trait, Model_Marker, Class, avgPA, -Model) %>%
    filter(!is.na(avgPA)) %>%  # Filter out NA values
    distinct() %>%
    spread(Model_Marker, avgPA)
  
  # Set rownames
  heatmap_data<-as.data.frame(heatmap_data)
  rownames(heatmap_data) <- heatmap_data$Trait
  heatmap_data <- heatmap_data[, -c(1, 2)]  # Remove Trait and Class columns
  # Convert to matrix
  heatmap_matrix <- as.matrix(heatmap_data)
  # Handle missing values
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  # Create color mappings
  unique_classes <- unique(merged_data$Class)
  class_colors <- brewer.pal(length(unique_classes), "Set3")
  color_mapping <- setNames(class_colors, unique_classes)
  
  # Apply color mapping
  function_colors <- class_data %>%
    select(Trait, Class) %>%
    distinct(Trait, Class) %>%
    mutate(color = color_mapping[Class])
  
  # Reorder rows in the heatmap matrix
  heatmap_matrix <- heatmap_matrix[match(function_colors$Trait, rownames(heatmap_matrix)), ]
  # Sort by the first column (e.g., "RR-CV_DMR")
  heatmap_matrix <- heatmap_matrix[order(heatmap_matrix[, "RR-REML_DMR"]), ]
  row_dendrogram <- hclust(dist(heatmap_matrix))
  # Reorder row_side_colors
  row_side_colors <- function_colors$color
  # Define the heatmap color palette
  heatmap_colors <- colorRampPalette(c("#9521f5", "#fbdbf5",  "#f40b56"))(100)
  
  # Generate and save the heatmap
  pdf(paste0(outPlot, "06.2_PA_heatmap_", class_name, ".pdf"))
  heatmap.2(heatmap_matrix,
            col = heatmap_colors,            # Color palette for heatmap
            trace = "none",                  # No trace lines
            dendrogram = "row",             # No dendrogram (manual sorting)
            RowSideColors = row_side_colors, # Add colors for each class
            margins = c(8, 8),               # Margins around the plot
            main = class_name,               # Title for the class
            xlab = "Marker", ylab = "Trait", # Axis labels
            key = TRUE, keysize = 1.2,       # Color key
            Rowv = as.dendrogram(row_dendrogram),                    # Disable row clustering
            Colv = FALSE,                    # Disable column clustering
            density.info = "none")           # No density plot
  dev.off()
}
########################
# compare PA between SNP and DMR per trait group
stat_results <- merged_data %>%
  filter(Model == "RR-REML") %>%
  group_by(Class) %>%
  summarise(p_value = t.test(avgPA[Marker == "SNP"], avgPA[Marker == "DMR"], paired = TRUE)$p.value)

# Print results
print(stat_results)
# Class                       p_value
# <chr>                         <dbl>
# 1 AA_and_NA_derivative       7.51e-10
# 2 Alkaloids_and_Terpenes     3.53e- 8
# 3 Flavonoids_and_derivatives 2.76e-13
# 4 Lipids                     9.67e- 3
# 5 Others                     4.07e- 9
# 6 Phenolamides               9.47e- 4
# 7 Polyphenols                1.94e- 4

# Merge p-values with filtered_data
filtered_data <- filtered_data %>%
  left_join(stat_results, by = "Class")

# Boxplot with p-values
ggplot(filtered_data, aes(x = Class, y = avgPA, fill = Marker)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Remove outliers for cleaner visualization
  scale_fill_manual(values = c("SNP" = "#999999", "DMR" = "#c06a80")) +
  labs(title = "PA Distribution by Trait Class",
       x = "Trait Class", y = "Predictive Ability (PA)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = stat_results, aes(x = Class, y = 1, 
                                     label = paste0("p = ", formatC(p_value, format = "e", digits = 2))),
            size = 4, vjust = 0, hjust = 0.5, inherit.aes = FALSE) 
ggsave(paste0(outPlot, "06.3_boxplot-trait.pdf"), width = 40, height = 20, units = "cm")

######################################################################
##### D Characterize Traits per H2 and PA
######################################################################
merged_data
binnedPA <- merged_data %>%
  mutate(PA_bin = cut(avgPA, 
                      breaks = seq(from = 0, # Rounds the min PA down to the nearest 0.1
                                   to = 1, # Set the upper limit as 1
                                   by = 0.1), # Increment by 0.1
                      include.lowest = TRUE))
binnedPA$association <- factor(binnedPA$association, levels = c("sig", "nonSig-GIF", "NonSig"))
mean_h2_data <- binnedPA %>%
  group_by(PA_bin, Marker) %>%
  summarise(mean_h2 = mean(h2, na.rm = TRUE)) %>%
  ungroup()

binnedPA
# Create the plot
ggplot() +
  geom_line(data = mean_h2_data, aes(x = PA_bin, y = mean_h2, color = Marker, 
                                     group = Marker), size = 2) +
  geom_point(data = binnedPA, aes(x = PA_bin, y = h2, color = Marker, 
                                  fill=Marker, alpha=0.1), 
             position = position_jitter(width = 0.1, height = 0)) +
  labs(title = "Mean h2 per PA Bin for Each Marker",
       x = "Predictive Ability (PA Range)",
       y = "Heritability (h2)",
       color = "Marker") +
  scale_color_manual(values=c("#c06a80", "#999999"))+
  scale_fill_manual(values=c("#c06a80","#999999"))+
  theme_minimal() +
  scale_x_discrete(guide = guide_axis(angle = 45)) +  # Rotate x-axis labels for readability
  facet_grid(Model ~ ., scales = "free")  # Separate plots by Model
ggsave(paste0(outPlot,"06.4_PA-H2.pdf"), width = 40, height = 20, units = "cm")
########
# Get paired SNP and DMR avgPA per Trait × Model
gain_table <- merged_data %>%
  select(Trait, Model, Marker, avgPA, h2) %>%
  pivot_wider(names_from = Marker, values_from = avgPA) %>%
  rename(PA_SNP = SNP, PA_DMR = DMR) %>%
  mutate(
    Relative_Gain = ifelse(abs(PA_SNP) > 1e-6,
                           ((PA_DMR - PA_SNP) / abs(PA_SNP)) * 100,
                           NA_real_),  # Avoid division by 0
    Fold_Change = ifelse(abs(PA_SNP) > 1e-6,
                         PA_DMR / PA_SNP,
                         NA_real_)
  )
# Overall across all traits and models
t_test_result <- t.test(
  x = gain_table$PA_DMR[!is.na(gain_table$PA_DMR) & !is.na(gain_table$PA_SNP)],
  y = gain_table$PA_SNP[!is.na(gain_table$PA_DMR) & !is.na(gain_table$PA_SNP)],
  paired = TRUE,
  alternative = "greater"
)
# For each model
model_stats <- gain_table %>%
  group_by(Model) %>%
  summarise(
    p_t = t.test(PA_DMR, PA_SNP, paired = TRUE, alternative = "greater")$p.value,
    p_wilcox = wilcox.test(PA_DMR, PA_SNP, paired = TRUE, alternative = "greater")$p.value,
    mean_gain = mean(Relative_Gain, na.rm = TRUE),
    median_gain = median(Relative_Gain, na.rm = TRUE),
    mean_fold = mean(Fold_Change, na.rm = TRUE),
    traits_better_with_DMR = sum(PA_DMR > PA_SNP, na.rm = TRUE),
    total_traits = n(),
    percent_better = mean(PA_DMR > PA_SNP, na.rm = TRUE) * 100
  )
# Filter only traits where all 3 models are present
DMR_better_traits <- gain_table %>%
  filter(!is.na(PA_DMR), !is.na(PA_SNP)) %>%
  mutate(DMR_better = PA_DMR > PA_SNP) %>%
  group_by(Trait) %>%
  summarise(DMR_better_all = all(DMR_better), .groups = "drop")%>%
  filter(DMR_better_all == TRUE)

# View or save
write.table(DMR_better_traits, file = paste0(outDir, "13.0_traits_DMR_better_all_models.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

snp_better_traits <- gain_table %>%
  filter(!is.na(PA_DMR), !is.na(PA_SNP)) %>%
  mutate(SNP_better = PA_SNP > PA_DMR) %>%
  group_by(Trait) %>%
  summarise(SNP_better_all = all(SNP_better), .groups = "drop") %>%
  filter(SNP_better_all == TRUE)

write.table(snp_better_traits, file = paste0(outDir, "13.1_traits_SNP_better_all_models.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


#
ggplot(model_stats, aes(x = reorder(Model, -mean_gain), y = mean_gain)) +
  geom_bar(stat = "identity", fill = "#c06a80", alpha = 0.8) +
  geom_text(aes(label = paste0("p=", formatC(p_t, format = "e", digits = 2))),
            vjust = -0.5, size = 4) +
  labs(title = "Mean Relative Gain in Predictive Ability (DMR vs SNP)",
       x = "Model",
       y = "Mean Relative Gain (%)") +
  theme_minimal() +
  ylim(0, max(model_stats$mean_gain, na.rm = TRUE) * 1.2)
ggsave(paste0(outPlot,"06.5_t-test.pdf"), width = 40, height = 20, units = "cm")

##
heatmap_data <- merged_data %>%
  ungroup() %>%  # <- THIS removes inherited grouping
  filter(Trait %in% DMR_better_traits$Trait) %>%
  mutate(Marker_Model = paste(Model, Marker, sep = "_")) %>%
  select(Trait, Marker_Model, avgPA) %>%
  distinct(Trait, Marker_Model, .keep_all = TRUE) %>%
  pivot_wider(names_from = Marker_Model, values_from = avgPA)

# Base R for rownames
rownames(heatmap_data) <- heatmap_data$Trait
heatmap_data <- heatmap_data[, -1]  # Remove Trait column

# Ensure matrix is clean
heatmap_matrix <- as.matrix(heatmap_data)
heatmap_matrix[is.na(heatmap_matrix)] <- 0  # replace NAs

# Try writing again
pdf(paste0(outPlot, "06.6_heatmap_DMR_better_traits.pdf"), width = 20, height = 14)
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("white", "#c06a80"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = 45,
         fontsize_row = 6,
         fontsize_col = 10,
         main = "DMR > SNP Traits (All Models)")
dev.off()
# Top 20 traits with highest DMR gain (RR-REML example)
top_traits <-gain_table %>%
  filter(Model == "RR-REML") %>%
  filter(!is.na(Relative_Gain),  # remove invalid
         !is.infinite(Relative_Gain),
         PA_SNP > 0.1,           # avoid division-by-zero bias
         PA_DMR > 0.5,           # focus on traits with good DMR performance
         Relative_Gain > 0) %>%
  ungroup() %>% 
  arrange(desc(Relative_Gain)) %>%
  slice_head(n = 20)
top_plot_data <- merged_data %>%
  filter(Trait %in% top_traits$Trait, Model == "RR-REML")

ggplot(top_plot_data, aes(x = Marker, y = avgPA, fill = Marker)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_wrap(~ Trait + Model, scales = "free_y") +
  scale_fill_manual(values = c("SNP" = "#999999", "DMR" = "#c06a80")) +
  labs(title = "Top 20 Traits with Highest Relative Gain from DMRs",
       y = "Predictive Ability", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.8, "lines"))

ggsave(paste0(outPlot, "06.7_top20_traits_DMR_gain.pdf"), width = 16, height = 10)


#
weighted_stats <- gain_table %>%
  filter(!is.na(h2)) %>%
  group_by(Model) %>%
  summarise(
    unweighted_gain = mean(Relative_Gain , na.rm = TRUE),
    weighted_gain = sum(Relative_Gain * h2, na.rm = TRUE) / sum(h2, na.rm = TRUE),
    weighted_fold = sum(Fold_Change * h2, na.rm = TRUE) / sum(h2, na.rm = TRUE)
  )



write.table(gain_table, file = paste0(outDir, "12.1_relative-gain_per_trait.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(model_stats, file = paste0(outDir, "12.2_model-level_gain_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(weighted_stats, file = paste0(outDir, "12.3_weighted_gain_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
