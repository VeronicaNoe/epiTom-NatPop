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
inDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/bc_predicted-values/"
outDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/results/"
outPlot <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/plots/"

########################## ##########################
# A Load and Process Data
####################################################
load_data <- function(file) {
  fread(file, sep = '\t', data.table = T, header = T, fill = TRUE, na.string = "NA", nThread = 20)
}
{
# Load Master Table
masterTable <- load_data(paste0(outDir, '08.1_master_table.tsv')) %>%
  filter(kinship == "DMR-SNP")

# traitGWAS <- load_data(paste0(outDir, '02.3_trait_H2.tsv')) %>%
#   filter(kin == "DMR-SNP") %>%
#   select(meta, h2, pseudo) %>%
#   distinct() %>%
#   arrange(meta)

# metaInfo <- load_data(paste0(outDir, '02.3_trait_H2.tsv')) %>%
#   filter(kin == "DMR-SNP") %>%
#   select(-h2, -pseudo, -kin) %>%
#   arrange(meta) %>%
#   group_by(meta) %>%
#   summarise(association = ifelse("sig" %in% association, "sig",
#                                  ifelse("nonSig-GIF" %in% association, "nonSig-GIF", "NonSig"))) %>%
#   left_join(traitGWAS, by = c("meta" = "meta")) %>%
#   arrange(meta)

#metaInfo <- metaName %>% left_join(metaInfo, by = c("Name" = "meta"))
}
########################################
# Load Prediction Data
########################################
#dataSet<-c("PCA-balance", "PCA-number-balance","random-number-balance","unbalance")
#PCA-number-balance: 80k-both
#PCA-balance: 120k-snp-80k-dmr.tsv
#random-number-balance: wo-na-1313-loci.tsv
#unbalance: snp-wo-na-all-dmrs.tsv  #original
metaName <- load_data("/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype/raw/metabolites_name.csv")
#
input <- list.files(path = paste0(inDir,"unbalance"), pattern = "LOCO", full.names = FALSE)
traitP <- do.call(rbind, lapply(input, function(i) load_data(paste0(inDir,"unbalance/", i))))
datasetType<-unlist(strsplit(input[1],'.', fixed = T))[1]
datasetType<-unlist(strsplit(datasetType,'_', fixed = T))[5]
#traitP$Cluster_Type<-NULL
traitP$Group<-"unbalance"
traitP$Dataset_Type<-datasetType
#

# Traits with avgPA < 0
toRemove<-traitP %>%
  group_by(Trait, Model, Marker, Dataset_Type) %>%
  filter(Model == "RR-REML") %>%
  summarise(all_na = all(is.na(PA)),
            .groups = "drop") %>%
  filter( all_na) %>%
  distinct(Trait) %>%
  pull(Trait)


# # Compute Mean Metrics
merged_data <- traitP %>%
  filter(!Trait %in% toRemove) %>%
  group_by(Trait, Model, Marker, Dataset_Type) %>%
  summarise(
    success_rate = mean(Success_Rate, na.rm = TRUE),
    avgPA = mean(PA, na.rm = TRUE),
    avgMSPE = mean(MSPE, na.rm = TRUE),
    avgMAPE = mean(MAPE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Trait, Model, Marker, Dataset_Type) %>%
  mutate(avgPA = ifelse(avgPA < 0, 0, avgPA)) %>%
  left_join(metaName, by = c("Trait" = "Name"))
#%>%
 #   left_join(metaInfo, by = c("Trait" = "Name"))
  #
#
pa_comparison <- merged_data %>%
  filter(Model == "RR-REML") %>%
  select(Trait,Class,`Compound Name`, Model,Marker, avgPA) %>%
  pivot_wider(names_from = Marker, values_from = avgPA) %>%
  mutate(DMR_better = DMR > SNP,
         PA_diff = DMR - SNP)
fwrite(pa_comparison, file = paste0(outDir,"14.0_pa_summary_rr-reml.tsv"), sep = "\t")

summary_stats <- pa_comparison %>%
  group_by(Model, Class) %>%
  #group_by(Model) %>%
  summarise(
    mean_DMR = mean(DMR, na.rm = TRUE),
    #sd_DMR = sd(DMR, na.rm = TRUE),
    mean_SNP = mean(SNP, na.rm = TRUE),
    #sd_SNP = sd(SNP, na.rm = TRUE),
    mean_DMR_SNP_diff = mean(PA_diff, na.rm = TRUE),
    #sd_diff = sd(PA_diff, na.rm = TRUE),
    FC_DMR_SNP = mean(DMR, na.rm=T)/mean(SNP, na.rm=T)
  )%>%
  arrange(desc(mean_DMR_SNP_diff))
# A tibble: 7 × 6
# Groups:   Model [1]
# Model   Class                   mean_DMR mean_SNP mean_DMR_SNP_diff FC_DMR_SNP
# <chr>   <chr>                      <dbl>    <dbl>             <dbl>      <dbl>
# 1 RR-REML AA_and_NA_derivative       0.323    0.140             0.183       2.30
# 2 RR-REML Lipids                     0.293    0.114             0.178       2.56
# 3 RR-REML Others                     0.322    0.174             0.147       1.84
# 4 RR-REML Alkaloids_and_Terpenes     0.297    0.173             0.124       1.72
# 5 RR-REML Polyphenols                0.323    0.206             0.117       1.57
# 6 RR-REML Phenolamides               0.252    0.143             0.109       1.76
# 7 RR-REML Flavonoids_and_derivat…    0.282    0.177             0.105       1.59

wilcox_test <- wilcox.test(pa_comparison$DMR, pa_comparison$SNP, paired = TRUE, alternative = "two.sided")
p_value <- wilcox_test$p.value

summary_stats %>%
  mutate(
    p_value = wilcox_test$p.value,
    statistic = wilcox_test$statistic
  )
# A tibble: 1 × 9
# Model   mean_DMR sd_DMR mean_SNP sd_SNP mean_diff sd_diff  p_value statistic
# <chr>      <dbl>  <dbl>    <dbl>  <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 RR-REML    0.300  0.161    0.166  0.155     0.134   0.145 2.48e-37     47348
  # Categorize H2
  # merged_data <- merged_data %>%
  #   mutate(H2_group = case_when(
  #     h2 <= 0.25 ~ "<= 0.25",
  #     h2 > 0.25 & h2 <= 0.75 ~ "> 0.25 & <= 0.75",
  #     h2 > 0.75 ~ "> 0.75"
  #   ))
###########################################################
#### B summarize models
##############################################################################
  # DMR summary with SD
  dmr_summary <- merged_data %>%
    filter(Marker == "DMR") %>%
    group_by(Model) %>%
    summarise(
      DMR = paste0(round(mean(avgPA, na.rm = TRUE), 3), " ± ", round(sd(avgPA, na.rm = TRUE), 3)),
      DMR_mean = mean(avgPA, na.rm = TRUE),
      DMR_sd = sd(avgPA, na.rm = TRUE),
      SR_DMR_mean = mean(success_rate, na.rm = TRUE) * 100,
      SR_DMR_sd = sd(success_rate, na.rm = TRUE) * 100
    )

  # SNP summary with SD
  snp_summary <- merged_data %>%
    filter(Marker == "SNP") %>%
    group_by(Model) %>%
    summarise(
      SNP = paste0(round(mean(avgPA, na.rm = TRUE), 3), " ± ", round(sd(avgPA, na.rm = TRUE), 3)),
      SNP_mean = mean(avgPA, na.rm = TRUE),
      SNP_sd = sd(avgPA, na.rm = TRUE),
      SR_SNP_mean = mean(success_rate, na.rm = TRUE) * 100,
      SR_SNP_sd = sd(success_rate, na.rm = TRUE) * 100
    )

  # Merge all results
  model_summary_full <- full_join(dmr_summary, snp_summary, by = "Model") %>%
    mutate(
      Success_Rate = rowMeans(cbind(SR_DMR_mean, SR_SNP_mean), na.rm = TRUE),
      Success_SD = apply(cbind(SR_DMR_sd, SR_SNP_sd), 1, function(x) sd(x, na.rm = TRUE))
    )

  print(model_summary_full, n = Inf)
  data.table::fwrite(model_summary_full, file = paste0(outDir,"13.0_model_summary_table.tsv"), sep = "\t")
  # RR-REML as the most stable and accurate
  # Despite testing several regression models, RR-REML consistently achieved both the highest predictive accuracy and the most robust success rate across markers, making it the most suitable choice for evaluating epigenomic predictive power.
  # “RR-REML stands out as the most reliable and performant model across both genetic and epigenetic data sources.”
  # A tibble: 6 × 4
  # Model   DMR           SNP           Success_Rate
  # <chr>   <chr>         <chr>                <dbl>
  # 1 ENET    0.203 ± 0.183 0.168 ± 0.16          67.5
  # 2 LASSO   0.192 ± 0.185 0.153 ± 0.157         63.5
  # 3 RR-CV   0.303 ± 0.164 0.223 ± 0.163         90.9
  # 4 RR-REML 0.3 ± 0.161   0.146 ± 0.143         97.9
  # 5 Ridge   0.303 ± 0.161 0.222 ± 0.163         91.2
  # 6 sENET   0.173 ± 0.172 0.153 ± 0.153         75.2
  # Combine DMR and SNP summaries
  summary_long <- bind_rows(
    dmr_summary %>%
      transmute(Model, Marker = "DMR",
                PA = DMR_mean,
                PA_sd = DMR_sd,
                Success_Rate = SR_DMR_mean,
                Success_SD = SR_DMR_sd),
    snp_summary %>%
      transmute(Model, Marker = "SNP",
                PA = SNP_mean,
                PA_sd = SNP_sd,
                Success_Rate = SR_SNP_mean,
                Success_SD = SR_SNP_sd)
  ) %>%
    mutate(
      ymin = pmax(Success_Rate - Success_SD, 0),
      ymax = pmin(Success_Rate + Success_SD, 100)
    )


  ggplot(summary_long, aes(x = reorder(Model, Success_Rate), y = Success_Rate, fill = Marker)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.2, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c("DMR" = "#c06a80", "SNP" = "#999999")) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = "Success Rate by Model and Marker",
         y = "Success Rate (%)", x = "Model", fill = "Marker Type") +
    theme_minimal(base_size = 14)
  ggsave(paste0(outPlot, "06.01_model-summary_success-rate_per-marker.pdf"),
         width = 40, height = 20, units = "cm")

  ggplot(summary_long, aes(x = reorder(Model, Success_Rate), y = PA, fill = Marker)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = PA - PA_sd, ymax = PA + PA_sd),
                  width = 0.2, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c("DMR" = "#c06a80", "SNP" = "#999999")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Model Performance by Marker Type",
         y = "Predictive Ability (Mean ± SD)",
         fill = "Marker",
         x = "Model") +
    theme_minimal(base_size = 14)
  ggsave(paste0(outPlot, "06.01_model-summary_PA_per-marker.pdf"),
         width = 40, height = 20, units = "cm")
######################################################################
##### C Characterize Traits per class using RR-REML
######################################################################
rr_reml_data <- merged_data %>%
  filter(Model == "RR-REML") %>%
  mutate(Class = factor(Class, levels = c(setdiff(unique(Class), "Others"), "Others")))
class_splits <- split(rr_reml_data, rr_reml_data$Class)
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
    ungroup() %>%
    select(Trait, Class) %>%
    distinct(Trait, Class) %>%
    mutate(color = color_mapping[Class])

  # Reorder rows in the heatmap matrix
  # Keep only those traits that are in the matrix
  function_colors <- function_colors %>%
    filter(Trait %in% rownames(heatmap_matrix))
  # Reorder matrix rows to match function_colors
  heatmap_matrix <- heatmap_matrix[match(function_colors$Trait, rownames(heatmap_matrix)), ]
  ##### this should work, maybe we shoul do this in bash
  #setdiff(function_colors$Trait, rownames(heatmap_matrix))
  #[1] "Tm147" "Tm191
  # Sort by the first column (e.g., "RR-CV_DMR")
  heatmap_matrix <- heatmap_matrix[order(heatmap_matrix[, "RR-REML_DMR"]), ]
  #row_dendrogram <- hclust(dist(heatmap_matrix))
  # Reorder row_side_colors
  row_side_colors <- function_colors$color
  # Define the heatmap color palette
  heatmap_colors <- colorRampPalette(c("#9521f5", "#fbdbf5",  "#f40b56"))(100)

  # Generate and save the heatmap
  # for future, the row dendrogram is not ok, low PA are cluster with higher PA
  pdf(paste0(outPlot, "06.02_PA_heatmap_", class_name, ".pdf"))
  heatmap.2(heatmap_matrix,
            col = heatmap_colors,            # Color palette for heatmap
            trace = "none",                  # No trace lines
            dendrogram = "none",             # No dendrogram (manual sorting)
            RowSideColors = row_side_colors, # Add colors for each class
            margins = c(8, 8),               # Margins around the plot
            main = class_name,               # Title for the class
            xlab = "Marker", ylab = "Trait", # Axis labels
            key = TRUE, keysize = 1.2,       # Color key
            Rowv = FALSE,                    # Disable row clustering
            Colv = FALSE,                    # Disable column clustering
            density.info = "none")           # No density plot
  dev.off()
}
######################################################################
##### D compare PA between SNP and DMR per trait group
######################################################################
stat_results <-rr_reml_data %>%
  group_by(Class) %>%
  summarise(p_value = t.test(avgPA[Marker == "SNP"], avgPA[Marker == "DMR"], paired = TRUE)$p.value) %>%
  arrange(p_value)
print(stat_results)
# A tibble: 7 × 2
#  Class                       p_value
#  <fct>                         <dbl>
#1 Others                     6.70e-16
#2 AA_and_NA_derivative       1.33e-15
#3 Flavonoids_and_derivatives 1.28e-12
#4 Alkaloids_and_Terpenes     2.13e- 9
#5 Phenolamides               4.78e- 8
#6 Polyphenols                6.12e- 7
#7 Lipids                     6.87e- 5
# Merge p-values with filtered_data
rr_reml_data <- rr_reml_data %>%
  left_join(stat_results, by = "Class")
# Boxplot with p-values
ggplot(rr_reml_data, aes(x = Class, y = avgPA, fill = Marker)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Remove outliers for cleaner visualization
  scale_fill_manual(values = c("SNP" = "#999999", "DMR" = "#c06a80")) +
  labs(title = "PA Distribution by Trait Class",
       x = "Trait Class", y = "Predictive Ability (PA)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = stat_results, aes(x = Class, y = 1,
                                     label = paste0("p = ", formatC(p_value, format = "e", digits = 2))),
            size = 4, vjust = 0, hjust = 0.5, inherit.aes = FALSE)
ggsave(paste0(outPlot, "06.03_boxplot-trait.pdf"), width = 40, height = 20, units = "cm")
######################################################################
##### E Compute Pearson correlation H2 and PA
######################################################################
cor_labels <- rr_reml_data %>%
  group_by(Model, Marker) %>%
  summarise(
    r = round(cor(h2, avgPA, method = "pearson", use = "complete.obs"), 2),
    .groups = "drop"
  ) %>%
  mutate(
    #label = paste0("r = ", r),
    label = ifelse(Marker == "DMR", paste0("PA-dmr|H2 r = ", r),
                   paste0("PA-snp|H2 r = ", r)),
    y = ifelse(Marker == "DMR", 0.95, 0.88)  # Stagger labels vertically
  )

ggplot(rr_reml_data, aes(x = h2, y = avgPA, color = Marker)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  geom_text(data = cor_labels, aes(x = 0.05, y = y, label = label, color = Marker),
            inherit.aes = FALSE, hjust = 0, size = 4) +
  scale_color_manual(values = c("DMR" = "#c06a80", "SNP" = "#999999")) +
  labs(title = "Predictive Ability vs Heritability",
       x = "Heritability (h²)",
       y = "Predictive Ability (PA)",
       color = "Marker") +
  theme_minimal(base_size = 14)
ggsave(paste0(outPlot,"06.04_PA-H2.pdf"), width = 40, height = 20, units = "cm")
######################################################################
##### F Calculate the gain power
######################################################################
# Get paired SNP and DMR avgPA per Trait × Model
gain_table <- merged_data %>%
  select(Trait,Class, `Compound Name`, Model, Marker, avgPA)%>%#, h2) %>%
  pivot_wider(names_from = Marker, values_from = avgPA) %>%
  rename(PA_SNP = SNP, PA_DMR = DMR) %>%
  mutate(
    Relative_Gain = ifelse(abs(PA_SNP) > 1e-6,
                           ((PA_DMR - PA_SNP) / abs(PA_SNP)) * 100,
                           NA_real_),  # Avoid division by 0
    Fold_Change = ifelse(abs(PA_SNP) > 1e-6,
                         PA_DMR / PA_SNP,
                         NA_real_)#,
    #Weighted_Gain = Relative_Gain * h2 / h2,
    #Weighted_Fold = Fold_Change * h2 / h2
  )
mean_gain_by_class <- gain_table %>%
  filter(Model=="RR-REML") %>%
  group_by(Class) %>%
  summarise(
    mean_relative_gain = mean(Relative_Gain, na.rm = TRUE),
    sd_relative_gain = sd(Relative_Gain, na.rm = TRUE),
    n = sum(!is.na(Relative_Gain)),
    .groups = "drop"
  )


write.table(gain_table, file = paste0(outDir, "12.1_relative-gain_per_trait.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Overall across all traits and models
t_test_result <- t.test(
  x = gain_table$PA_DMR[!is.na(gain_table$PA_DMR) & !is.na(gain_table$PA_SNP)],
  y = gain_table$PA_SNP[!is.na(gain_table$PA_DMR) & !is.na(gain_table$PA_SNP)],
  paired = TRUE,
  alternative = "greater"
)

# For each model
# First summarize as before
trait_comparison_by_model <- gain_table %>%
  mutate(comparable = ifelse(!is.na(PA_DMR) & !is.na(PA_SNP) & PA_SNP > 0 & PA_DMR > 0, TRUE, FALSE)) %>%
  group_by(Model, comparable) %>%
  summarise(n_traits = n_distinct(Trait), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = comparable,
    values_from = n_traits,
    names_prefix = ""
  ) %>%
  rename(
    comparable = `TRUE`,
    non_comparable = `FALSE`
  ) %>%
  mutate(
    comparable = replace_na(comparable, 0),
    non_comparable = replace_na(non_comparable, 0)
  )
trait_comparison_by_model
# A tibble: 6 × 3
# Model   non_comparable comparable
# <chr>            <int>      <int>
# 1 ENET               156        183
# 2 LASSO              169        170
# 3 RR-CV               62        277
# 4 RR-REML             85        254
# 5 Ridge               59        280
# 6 sENET              153        186

model_stats <- gain_table %>%
  filter(!is.na(PA_DMR), !is.na(PA_SNP),
         PA_SNP > 0, PA_DMR > 0) %>%  # Optional: exclude zeros
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

weighted_stats <- gain_table %>%
  filter(!is.na(h2), !is.na(PA_DMR), !is.na(PA_SNP),
         PA_SNP > 0, PA_DMR > 0) %>%
  group_by(Model) %>%
  summarise(
    unweighted_gain = mean(Relative_Gain , na.rm = TRUE),
    weighted_gain = sum(Relative_Gain * h2, na.rm = TRUE) / sum(h2, na.rm = TRUE),
    weighted_fold = sum(Fold_Change * h2, na.rm = TRUE) / sum(h2, na.rm = TRUE)
  )

# Combine the two tables into one:
full_model_summary <- model_stats %>%
  left_join(weighted_stats, by = "Model")
write.table(full_model_summary, file = paste0(outDir, "12.2_model-level_gain_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(full_model_summary, n=Inf, width=Inf)
# A tibble: 6 × 12
# Model        p_t p_wilcox mean_gain median_gain mean_fold
# <chr>      <dbl>    <dbl>     <dbl>       <dbl>     <dbl>
#   1 ENET    1.33e- 2 1.75e- 2      235.        11.3      3.35
# 2 LASSO   6.45e- 3 7.67e- 3      351.        22.1      4.51
# 3 RR-CV   4.55e-21 2.42e-19      239.        34.1      3.39
# 4 RR-REML 5.43e-41 4.52e-32      569.        95.7      6.69
# 5 Ridge   1.87e-19 6.89e-18      184.        32.3      2.84
# 6 sENET   8.66e- 2 7.06e- 2      104.        10.8      2.04
# traits_better_with_DMR total_traits percent_better unweighted_gain
# <int>        <int>          <dbl>           <dbl>
#   1                    104          183           56.8            235.
# 2                    105          170           61.8            351.
# 3                    204          277           73.6            239.
# 4                    212          254           83.5            569.
# 5                    203          280           72.5            184.
# 6                    103          186           55.4            104.
# weighted_gain weighted_fold
# <dbl>         <dbl>
#   1         204.           3.04
# 2         259.           3.59
# 3         190.           2.90
# 4         626.           7.26
# 5         170.           2.70
# 6          92.4          1.92

# Filter only traits where all models are present
DMR_better_traits <- gain_table %>%
  filter(!is.na(PA_DMR), !is.na(PA_SNP)) %>%
  mutate(DMR_better = PA_DMR > PA_SNP) %>%
  group_by(Trait) %>%
  summarise(DMR_better_all = all(DMR_better), .groups = "drop")%>%
  filter(DMR_better_all == TRUE)
# 105
# View or save
write.table(DMR_better_traits, file = paste0(outDir, "13.0_traits_DMR_better_all_models.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

snp_better_traits <- gain_table %>%
  filter(!is.na(PA_DMR), !is.na(PA_SNP)) %>%
  mutate(SNP_better = PA_SNP > PA_DMR) %>%
  group_by(Trait) %>%
  summarise(SNP_better_all = all(SNP_better), .groups = "drop") %>%
  filter(SNP_better_all == TRUE)
# 19
write.table(snp_better_traits, file = paste0(outDir, "13.1_traits_SNP_better_all_models.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
## boxplot
gain_table_plot <- gain_table %>%
  filter(!is.na(Relative_Gain), is.finite(Relative_Gain)) %>%
  left_join(model_stats %>% select(Model, p_t), by = "Model")

ggplot(gain_table_plot, aes(x = reorder(Model, Relative_Gain, median, na.rm = TRUE),
                            y = Relative_Gain, fill = Model)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
  geom_text(data = model_stats,
            aes(x = Model, y = max(gain_table_plot$Relative_Gain, na.rm = TRUE) * 1.05,
                label = paste0("p = ", formatC(p_t, format = "e", digits = 2))),
            inherit.aes = FALSE, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Distribution of Relative Gain in Predictive Ability (DMR vs SNP)",
       y = "Relative Gain (%)",
       x = "Model") +
  theme_minimal(base_size = 14) +
  coord_cartesian(
    ylim = quantile(gain_table_plot$Relative_Gain, c(0.05, 0.95), na.rm = TRUE)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")+
  theme(legend.position = "none")
ggsave(paste0(outPlot,"06.05_power-gain_t-test_per-model_boxplot.pdf"), width = 40, height = 20, units = "cm")

##
#I don't know if this is informativa
{
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
pdf(paste0(outPlot, "06.06_heatmap_DMR_better_traits.pdf"), width = 20, height = 14)
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("white", "#c06a80"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = 45,
         fontsize_row = 6,
         fontsize_col = 10,
         main = "DMR > SNP Traits (All Models)")
invisible(dev.off())
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
ggsave(paste0(outPlot, "06.07_top20_traits_DMR_gain_RR-REML.pdf"), width = 16, height = 10)
}
##################
# Initialize results list
##################
  enrichment_list <- list()
  # Unique trait classes to test
  for (i in unique(merged_data$Class)) {
    # All traits tested in RR-REML with valid PA
    valid_traits <- gain_table %>%
      filter(Model == "RR-REML", !is.na(PA_DMR), !is.na(PA_SNP),
             PA_SNP > 0 & PA_DMR > 0) %>%
      distinct(Trait) %>%
      pull(Trait)

    # Get traits in the current class
    class_traits <- merged_data %>%
      filter(Model == "RR-REML", Trait %in% valid_traits,
             Marker == "DMR",Class == i) %>%
      distinct(Trait) %>%
      pull(Trait)

    # a: DMR-better & in class
    a <- gain_table %>%
      filter(Model == "RR-REML", Trait %in% class_traits,
             PA_DMR > PA_SNP, PA_SNP > 0 & PA_DMR > 0) %>%
      distinct(Trait) %>%
      nrow()

    # b: Not DMR-better & in class
    b <-  gain_table %>%
      filter(Model == "RR-REML", Trait %in% class_traits,
             PA_DMR < PA_SNP, PA_SNP > 0 & PA_DMR > 0) %>%
      distinct(Trait) %>%
      nrow()

    # c: DMR-better & NOT in class
    c <- gain_table %>%
      filter(Model == "RR-REML", Trait %in% setdiff(valid_traits, class_traits),
             PA_DMR > PA_SNP, PA_SNP > 0 & PA_DMR > 0) %>%
      distinct(Trait) %>%
      nrow()

    # d: Not DMR-better & NOT in class
    d <- gain_table %>%
      filter(Model == "RR-REML", Trait %in% setdiff(valid_traits, class_traits), PA_DMR <= PA_SNP, PA_SNP > 0 & PA_DMR > 0) %>%
      distinct(Trait) %>%
      nrow()

    # Fisher test
    contingency <- matrix(c(a,c,b, d), nrow = 2)
    pval <- fisher.test(contingency)$p.value

    # Store result
    enrichment_list[[i]] <- data.frame(
      Class = i,
      a = a, b = b, c = c, d = d,
      p_value = pval,
      prop_in_class = a / (a + b)
    )
  }
  # Combine and adjust p-values
  enrichment_results <- do.call(rbind, enrichment_list) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
    arrange(p_adj)
  enrichment_results
  # Class  a  b   c  d     p_value
  # AA_and_NA_derivative             AA_and_NA_derivative 37  0 175 42 0.001249769
  # Flavonoids_and_derivatives Flavonoids_and_derivatives 52 20 160 22 0.004393719
  # Lipids                                         Lipids 11  4 201 38 0.283270718
  # Others                                         Others 43  7 169 35 0.675783596
  # Alkaloids_and_Terpenes         Alkaloids_and_Terpenes 24  6 188 36 0.601861698
  # Polyphenols                               Polyphenols 26  3 186 39 0.434368940
  # Phenolamides                             Phenolamides 19  2 193 40 0.542975498
  # prop_in_class       p_adj
  # AA_and_NA_derivative           1.0000000 0.008748384
  # Flavonoids_and_derivatives     0.7222222 0.015378015
  # Lipids                         0.7333333 0.660965009
  # Others                         0.8600000 0.675783596
  # Alkaloids_and_Terpenes         0.8000000 0.675783596
  # Polyphenols                    0.8965517 0.675783596
  # Phenolamides                   0.9047619 0.675783596

  ## Result
  # We found a significant enrichment of traits from the Flavonoids and
  # derivatives class among those with higher predictive power using DMRs
  # compared to SNPs (Fisher’s exact test, FDR-adjusted p = 0.027).
  # Specifically, 66 of 92 flavonoid-related traits (~71.7%) showed higher
  # predictive ability when modeled using DMRs. This suggests that DMRs may
  # better capture epigenetic variation relevant to flavonoid metabolism.
####################################
# Prepare gain table with class info
####################################

# STEP 1: Join gain data with class info (only RR-REML, valid comparisons)
  # Total traits per class (based on RR-REML + DMR entries)
  valid_traits <- gain_table %>%
    filter(Model == "RR-REML", !is.na(PA_DMR), !is.na(PA_SNP),
           PA_SNP > 0 & PA_DMR > 0) %>%
    distinct(Trait) %>%
    pull(Trait)

  gain_with_class <- gain_table %>%
    ungroup() %>%
    filter(Model == "RR-REML", Trait %in% valid_traits,
           PA_DMR > PA_SNP, PA_SNP > 0 & PA_DMR > 0,
           #!is.na(h2),
           !is.na(PA_DMR), !is.na(PA_SNP),
           !is.na(Relative_Gain), is.finite(Relative_Gain))
  #212
  # Get total number of traits per class (for % calculation)
  total_traits_by_class <- merged_data %>%
    ungroup() %>%
    filter(Model == "RR-REML", Marker == "DMR",Trait %in% valid_traits) %>%
    select(Trait, Class) %>%
    distinct() %>%
    count(Class, name = "total_traits")

  # Calculate summary stats per class
  class_bubble_data <- gain_with_class %>%
    filter(!is.na(Relative_Gain), is.finite(Relative_Gain)) %>%
    group_by(Class) %>%
    summarise(
      n_traits = n(),  # number of traits better predicted by DMR
      median_relative_gain = median(Relative_Gain, na.rm = TRUE),
      mean_fold_change = mean(Fold_Change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(total_traits_by_class, by = "Class") %>%
    mutate(
      percent_dmr_better = round((n_traits / total_traits) * 100, 1)
    )%>%
    mutate(
      fold_bin = cut(mean_fold_change,
                     breaks = c(0, 2, 3, 4, 5, 10, Inf),
                     labels = c("≤2", "2", "3", "4", "5", ">10"),
                     include.lowest = TRUE),
      trait_bin = cut(n_traits,
                      breaks = c(0, 20, 40, Inf),
                      labels = c("Small", "Medium", "Large"),
                      include.lowest = TRUE)
    )%>%
    arrange(desc(mean_fold_change))
  class_bubble_data
  # Bubble plot
  ggplot(class_bubble_data, aes(
    x = median_relative_gain,
    y = percent_dmr_better,
    size = trait_bin,
    fill = fold_bin
  )) +
    geom_point(shape = 21, color = "black", alpha = 0.9) +
    geom_text(aes(label = Class), vjust = -1, size = 4) +
    scale_fill_brewer(palette = "YlOrRd", name = "Mean Fold\nChange (binned)") +
    scale_size_manual(values = c(Small = 5, Medium = 8, Large = 11),
                      name = "# Traits") +
    labs(
      title = "DMR Prediction Advantage by Trait Class",
      x = "Median Relative Gain (%)",
      y = "% of Traits Where DMR > SNP"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right")
  # Save
  ggsave(paste0(outPlot, "06.08_summary_class_bubbleplot.pdf"), width = 12, height = 8)
  # Among all trait classes tested using RR-REML, flavonoids showed the strongest enrichment
  # for DMR-based predictive improvement, with over 79% of traits having higher PA using DMRs.
  # The median relative gain across trait classes ranged from ~30% to 90%, with a typical fold-change
  # of about 2 to 4.
