suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(data.table)
  library(gprofiler2)
  library(tidyr)
  library(ggpattern)
  library(tidyverse)
  library(purrr)
})
setwd("/Users/macbook/Desktop/tmo/Fig4")

# setwd("/mnt/disk2/vibanez/10_data-analysis/Fig4/results")
# outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/results/"
# outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/plots/"

################################################################################
## combine GE for genes with change in the teM frequency and KOs
################################################################################
# function
clean_revigo <- function(file_path) {
  revigo <- data.table::fread(file_path,
                              sep = "\t",
                              data.table = TRUE,
                              fill = TRUE,
                              check.names = FALSE,
                              na.strings = c("NA", ""))
  
  new_term_ids <- character(nrow(revigo))
  new_term_names <- character(nrow(revigo))
  current_rep <- NA_character_
  current_name <- NA_character_
  
  for (i in seq_len(nrow(revigo))) {
    if (revigo$Representative[i] != "null") {
      new_term_ids[i] <- current_rep
      new_term_names[i] <- current_name
    } else {
      new_term_ids[i] <- revigo$TermID[i]
      new_term_names[i] <- revigo$Name[i]
      current_rep <- revigo$TermID[i]
      current_name <- revigo$Name[i]
    }
  }
  
  revigo$new_term_id <- new_term_ids
  revigo$new_term_name <- new_term_names
  
  revigo_clean <- revigo %>%
    dplyr::select(term_id = TermID, new_term_id, new_term_name)
  
  return(revigo_clean)
}

#
DEG_genes_teM_freq_change <- data.table::fread(
  "05.05.DEG_Pim-SLC_teM-frequency-change.tsv",#paste0(outDir, "07.01.GO_freq-change-PIM-SLC_DEG.tsv"),
  sep = ",",
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = c("NA", "")  
)
teM_freq_change <- data.table::fread(
  "03.02.accGroups_teM-frequency_PIM-SLC.tsv",#paste0(outDir, "07.01.GO_freq-change-PIM-SLC_DEG.tsv"),
  sep = "\t",
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = c("NA", "")  
)
bg_genes<-unique(teM_freq_change$geneName)
# get sample Gene expression
comparison<-c("SLC", "PIM","DEG")
all_go_results <- list()  # Initial
for(i in comparison){
  if(i=="DEG"){
    top_genes_long <- DEG_genes_teM_freq_change %>%
      filter(te_significance=="downregulated" | !is.na(pe_significant_direction)) %>%
      pull(geneName)%>%
      unique()
    source<-"DEG_genes_teM_freq_change"
  }else{
    top_genes_long <- teM_freq_change %>%
      filter(higher_in==i, sig_color!="grey") %>%
      pull(geneName)%>%
      unique()
    source<-paste0(i,"_higher")
  }
  # Run g:Profiler query
  res <- gost(query = top_genes_long, 
              organism = "slycopersicum",  # Code for Solanum lycopersicum
              sources = c("GO:BP"),
              domain_scope = "custom",
              evcodes = TRUE,
              custom_bg = bg_genes,
              significant = FALSE)
 
  # Access the table
  go_results <- res$result
  go_results$source<-source
  go_results<-go_results[, c("term_id", "term_name", "intersection_size", "term_size","intersection","effective_domain_size","source")]
  all_go_results[[i]] <- go_results
}

all_go_results <- rbindlist(all_go_results)
head(all_go_results)
length(unique(all_go_results$term_name))
#res$meta$version
#"e112_eg59_p19_25aa4782"
###
# Calculate proportion first
all_go_results <- all_go_results %>%
  mutate(percentage = (intersection_size*100 / term_size )) 

# Pivot wider for intersection_size
intsize_wide <- all_go_results %>%
  select(term_id, term_name, term_size,source, intersection_size) %>%
  pivot_wider(names_from = source, values_from = intersection_size, names_prefix = "int_", values_fill = 0)%>%
  filter(int_DEG_genes_teM_freq_change!=0)

data.table::fwrite(intsize_wide, file = "07.01.GO_freq-change-PIM-SLC_DEG.tsv") #clean annoaiton

#
genes_wide <- all_go_results %>%
  select(term_id, term_name, source, intersection) %>%
  pivot_wider(names_from = source, values_from = intersection, names_prefix = "genes_", values_fill =  list(intersection = NA_character_)) %>%
  filter(genes_DEG_genes_teM_freq_change!=0)

# Combine both
go_combined <- left_join(intsize_wide, genes_wide, by = c("term_id", "term_name"))

go_combined

###########################
# http://revigo.irb.hr/
# default + small
revigo_cleaned <- clean_revigo("Revigo_BP_Table.tsv")

go_long<-all_go_results %>%
  filter(term_id %in% go_combined$term_id) %>%
  select(term_id, term_name, source, intersection, intersection_size, term_size) %>%
  left_join(revigo_cleaned, by = "term_id") %>%
  mutate(merged_id = ifelse(is.na(new_term_id), term_id, new_term_id),
         merged_name = ifelse(is.na(new_term_name), term_name, new_term_name))

go_long_separated <- go_long %>%
  separate_rows(intersection, sep = ",") %>%
  filter(intersection != "")
go_grouped <- go_long_separated %>%
  group_by(merged_id, merged_name, source) %>%
  summarise(
    merged_term_size = max(term_size),  # conservative estimate
    merged_intersection_size = n_distinct(intersection),
    merged_genes = paste(unique(intersection), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(percentage = (merged_intersection_size*100) / merged_term_size)

go_wide_final <- go_grouped %>%
  select(merged_id, merged_name, merged_term_size,source, percentage, merged_intersection_size, merged_genes) %>%
  pivot_wider(
    names_from = source,
    values_from = c(percentage, merged_intersection_size, merged_genes),
    values_fill = list(percentage = 0, merged_intersection_size = 0, merged_genes = "")
  ) %>%
  arrange(desc(merged_term_size))%>%
  filter(merged_term_size < 240)


# Create helper function to split gene lists
split_genes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  str_split(x, ",")[[1]] %>% str_trim()
}

# Expand per GO term
plot_data <- go_wide_final %>%
  rowwise() %>%
  mutate(
    genes_DEG = list(split_genes(merged_genes_DEG_genes_teM_freq_change)),
    genes_PIM = list(split_genes(merged_genes_PIM_higher)),
    genes_SLC = list(split_genes(merged_genes_SLC_higher)),
    
    deg_pim = length(intersect(genes_DEG, genes_PIM)),
    deg_slc = length(intersect(genes_DEG, genes_SLC)),
    
    pim_higher = length(genes_PIM),
    slc_higher = length(genes_SLC)
  ) %>%
  ungroup() %>%
  select(merged_name, merged_term_size, deg_pim, deg_slc, pim_higher, slc_higher) %>%
  pivot_longer(cols = -c(merged_name, merged_term_size), names_to = "group", values_to = "count") %>%
  mutate(percentage = 100 * count / merged_term_size)

# Make sure pattern variable exists
plot_data <- plot_data %>%
  mutate(group_type = case_when(
    str_detect(group, "deg_") ~ "DEG",
    TRUE ~ "Higher"
  ))

# Order by SLC_higher
ordering <- plot_data %>%
  filter(group == "slc_higher") %>%
  arrange(desc(percentage)) %>%
  pull(merged_name)

# Plot with ggpattern
ggplot(plot_data, aes(x = factor(merged_name, levels = ordering), y = percentage,
                      fill = group, pattern = group_type)) +
  geom_bar_pattern(
    stat = "identity", position = position_dodge(width = 0.8),
    color = "black",
    pattern_density = 0.4,
    pattern_spacing = 0.03,
    pattern_fill = "black",
    pattern_angle = 45
  ) +
  scale_fill_manual(values = c(
    "pim_higher" = "#00a681ff",
    "deg_pim"    = "#00a681ff",
    "slc_higher" = "#ed83b5ff",
    "deg_slc"    = "#ed83b5ff"
  )) +
  scale_pattern_manual(values = c("Higher" = "none", "DEG" = "stripe")) +
  labs(x = "GO Term", y = "% of genes", fill = "Group", pattern = "Type") +
  coord_flip() +
  theme_minimal(base_size = 12)

ggsave("07.02_shared_GO_terms.pdf", width = 20, height = 20, units = "cm")


########################
### Get expression per GO term
########################
sampleGE<-data.table::fread(
  "05.01.epiallele-geneExpression.tsv", #paste0(outDir,"05.01.epiallele-geneExpression.tsv"), 
  sep = ',', data.table = T, fill = TRUE, check.names=FALSE,
  na.string=c("NA"), nThread = 10)

# get KO gene expression
#inDir<-"/mnt/disk2/vibanez/09_KO-processing/09.1_rnaseq/ac_counts"
KOsGE<-data.table::fread("00.all-KOs_counts_log2FC.tsv", sep = '\t', 
                         data.table = T, fill = TRUE, check.names=FALSE,
                         na.string=c("NA"), nThread = 10)

###
# Combine both tables
desired_order <- c(
  "PIM_UM",  "PIM_teM",
  "SLC_UM",  "SLC_teM", 
  "SLL_UM",  "SLL_teM", 
  "ddm1", "wt-ddm1"
)

DEG_go<-go_summary_top_clean %>%
  filter(source == "DEG_genes_teM_freq_change")

# 1. Get all genes from DEG terms
gene_list <- DEG_go$intersection %>%
  str_split(",") %>%
  unlist() %>%
  unique()

# 2. Filter expression data
filtered_sampleGE <- sampleGE %>%
  filter(geneName %in% gene_list, sampleEpiallele != 'gbM', accGroup != "WILD") %>%
  mutate(group = paste(accGroup, sampleEpiallele, sep = "_"))

filtered_KOsGE <- KOsGE %>%
  filter(Geneid %in% gene_list, KO == "ddm1")

# 3. Prepare data for plotting
sampleGE_clean <- filtered_sampleGE %>%
  transmute(geneName, expressionLog2, group)

KOsGE_clean <- filtered_KOsGE %>%
  transmute(geneName = Geneid, expressionLog2 = expLog2, group = sample)

combined_GE <- bind_rows(sampleGE_clean, KOsGE_clean) %>%
  mutate(group = factor(group, levels = desired_order)) %>%
  mutate(epiallele_group = case_when(
    str_detect(group, "teM") ~ "teM",
    str_detect(group, "UM") ~ "UM",
    group == "ddm1" ~ "UM",
    group == "wt-ddm1" ~ "wt"
  ))

epiallele_colors <- c("teM" = "#820a86", "UM" = "grey", "wt" = "black")

# 4. Plot all genes together (faceted)
count_table <- combined_GE %>%
  group_by(geneName, group) %>%
  summarise(n = n(), .groups = "drop")

go_terms_by_gene <- DEG_go %>%
  separate_rows(intersection, sep = ",") %>%
  mutate(intersection = str_trim(intersection)) %>%
  filter(intersection %in% gene_list) %>%
  select(geneName = intersection, term_id, term_name, source) %>%
  group_by(geneName) %>%
  summarise(go_terms = paste(sort(unique(term_name)), collapse = "\n"), .groups = "drop")

combined_GE_annot <- combined_GE %>%
  left_join(go_terms_by_gene, by = "geneName") %>%
  mutate(panel_label = paste0(geneName, "\n", go_terms))  # New label

ggplot(combined_GE_annot, aes(x = group, y = expressionLog2, fill = epiallele_group)) +
  geom_jitter(aes(color = epiallele_group), width = 0.2, size = 1.2, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.7) +
  labs(title = paste0("Expression of DEG genes (", length(gene_list), " total)"),
       x = "Group", y = "Expression (log2)") +
  scale_fill_manual(values = epiallele_colors) +
  scale_color_manual(values = epiallele_colors) +
  facet_wrap(~ panel_label, scales = "free_y") +  # Use GO-annotated labels
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10),
        legend.position = "none")
ggsave("06.05.DEG_set_combined_plot_with_GO-term.pdf", width = 20, height = 20, units = "cm")

ggplot(combined_GE_annot, aes(x = group, y = expressionLog2, fill = epiallele_group)) +
  geom_jitter(aes(color = epiallele_group), width = 0.2, size = 1.2, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.7) +
  geom_text(data = count_table, aes(x = group, y = max(combined_GE$expressionLog2, na.rm = TRUE) + 0.5, label = n),
            inherit.aes = FALSE, size = 3, color = "black") +
  labs(title = paste0("Expression of DEG genes (", length(gene_list), " total)"),
       x = "Group", y = "Expression (log2)") +
  scale_fill_manual(values = epiallele_colors) +
  scale_color_manual(values = epiallele_colors) +
  facet_wrap(~ geneName, scales = "free_y") +  # Use GO-annotated labels
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10),
        legend.position = "none")
ggsave("06.05.DEG_set_combined_plot.pdf", width = 20, height = 20, units = "cm")


# Flatten and get unique gene names
gene_list <- unique(unlist(DEG_go$all_genes))
# Convert to data.table and write to file
data.table::fwrite(
  data.table::data.table(gene_name = gene_list),
  file = "gene_list.tsv",
  sep = "\t",
  quote = FALSE, col.names = F
)
# get methylation over tose genes
######################
###
######################
# Set directories
# Save result tables
all_samples_df<-data.table::fread("06.07.all_methylation_levels_long.tsv", 
                            sep = '\t', data.table = T, fill = TRUE, check.names=FALSE,
                            na.string=c("NA"), nThread = 10)

# get samples and it's expression
sampleGE<-data.table::fread("05.01.epiallele-geneExpression.tsv", 
                            sep = ',', data.table = T, fill = TRUE, check.names=FALSE,
                            na.string=c("NA"), nThread = 10)

sampleGE_clean <- sampleGE %>%
  filter(sampleEpiallele!="gbM", accGroup!="WILD") %>%
  select(geneName, accessions, sampleEpiallele, accGroup) %>%
  rename(
    gene = geneName,
    sample = accessions
  )

samples2keep<-intersect(sampleGE_clean$sample, all_samples_df$sample)
koSamples<-grep('ddm1', unique(all_samples_df$sample), value = T)
samples2keep<-c(samples2keep, koSamples)
# Join with all_samples_df
desired_order <- c(
  "PIM_UM",  "PIM_teM",
  "SLC_UM",  "SLC_teM", 
  "SLL_UM",  "SLL_teM",
  "ddm1", "wt-ddm1"
)

all_samples_annotated <- all_samples_df %>%
  filter(sample %in% samples2keep) %>%
  left_join(sampleGE_clean, by = c("gene", "sample")) %>%
  mutate(
    clean_sample = if_else(
      str_detect(sample, "^ddm1|^wt-ddm1"),
      str_remove(sample, "-[0-9]+$"),
      sample
    ),
    group = if_else(
      is.na(accGroup) | is.na(sampleEpiallele),
      clean_sample,
      paste(accGroup, sampleEpiallele, sep = "_")
    ),
    group = factor(group, levels = desired_order),
    epiallele_group = case_when(
      str_detect(group, "teM") ~ "teM",
      str_detect(group, "UM") ~ "UM",
      group == "ddm1" ~ "UM",       # ← this is what you meant
      group == "wt-ddm1" ~ "wt",
      TRUE ~ "other"
    )
  )
all_samples_annotated <- all_samples_annotated %>%
  left_join(go_terms_by_gene, by = c("gene" = "geneName")) %>%
  mutate(panel_label = paste0(gene, "\n", str_wrap(go_terms, width = 50)))

epiallele_colors <- c(
  "teM" = "#820a86",
  "UM"  = "grey",
  "wt" = "black"
)

count_table <- all_samples_annotated %>%
  group_by(gene, group) %>%
  summarise(n = n(), y_max = max(methylation_level, na.rm = TRUE), .groups = "drop")%>%
  mutate(label_y = y_max + 0.05)  # Shift text above max

ggplot(all_samples_annotated, aes(x = group, y = methylation_level, fill = epiallele_group)) +
  geom_jitter(aes(color = epiallele_group), width = 0.2, size = 1.2, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.7) +
  geom_text(data = count_table, aes(x = group, y = label_y, label = n),
            inherit.aes = FALSE, size = 3, color = "black") +
  labs(
    title = paste0("MethLevels of DEG genes (", length(unique(all_samples_annotated$gene)), " total)"),
    x = "Group", y = "Methylation level"
  ) +
  scale_fill_manual(values = epiallele_colors) +
  scale_color_manual(values = epiallele_colors) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    panel.spacing = unit(1.5, "lines"),
    legend.position = "none"
  )
ggsave("06.06.methLevel_selected-genes.pdf", width = 20, height = 20, units = "cm")

ggplot(all_samples_annotated, aes(x = group, y = methylation_level, fill = epiallele_group)) +
  geom_jitter(aes(color = epiallele_group), width = 0.2, size = 1.2, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.7) +
  #geom_text(data = count_table, aes(x = group, y = label_y, label = n),
  #          inherit.aes = FALSE, size = 3, color = "black") +
  labs(
    title = paste0("MethLevels of DEG genes (", length(unique(all_samples_annotated$gene)), " total)"),
    x = "Group", y = "Methylation level"
  ) +
  scale_fill_manual(values = epiallele_colors) +
  scale_color_manual(values = epiallele_colors) +
  facet_wrap(~ panel_label, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    panel.spacing = unit(1.5, "lines"),
    legend.position = "none"
  )
ggsave("06.06.methLevel_selected-genes_GO-term.pdf", width = 20, height = 20, units = "cm")



