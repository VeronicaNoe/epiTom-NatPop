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

gprofiler2::set_base_url("https://biit.cs.ut.ee/gprofiler")

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
PE_genes <- data.table::fread(
  "02.07_epiallele_sample.tsv",#paste0(outDir, "07.01.GO_freq-change-PIM-SLC_DEG.tsv"),
  sep = "\t",
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = c("NA", "")  
)

PE_genes<-PE_genes %>%
  filter(geneEpiallele == "PE") %>%
  select(geneName) %>%
  distinct(geneName) %>%
  pull(geneName)
bg_genes <- data.table::fread(
  "allGenes.bed",#paste0(outDir, "07.01.GO_freq-change-PIM-SLC_DEG.tsv"),
  sep = "\t",
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = c("NA", "")  
)
bg_genes<-bg_genes %>%
  filter(V1 != "SL2.50ch00") %>%
  select(V5) %>%
  distinct(V5) %>%
  pull(V5)

# Run g:Profiler query
res <- gost(query = PE_genes, 
              organism = "slycopersicum",  # Code for Solanum lycopersicum
              sources = c("GO:BP"),
              domain_scope = "custom",
              evcodes = TRUE,
              custom_bg = bg_genes,
              significant = FALSE)
 
# Access the table
go_results <- res$result
go_results$source<-"pe_genes"
go_results<-go_results[, c("term_id", "term_name", "intersection_size", "term_size","intersection","effective_domain_size","source")]
  
#res$meta$version
#"e112_eg59_p19_25aa4782"
###
# Calculate proportion first
all_go_results <- go_results %>%
  mutate(percentage = (intersection_size*100 / term_size )) 

# Pivot wider for intersection_size
intsize_wide <- all_go_results %>%
  select(term_id, term_name, term_size,source, intersection_size) %>%
  pivot_wider(names_from = source, values_from = intersection_size, names_prefix = "int_", values_fill = 0)

data.table::fwrite(intsize_wide, file = "07.01.GO_PE-genes.tsv") #clean annoaiton

#
genes_wide <- all_go_results %>%
  select(term_id, term_name, source, intersection) %>%
  pivot_wider(names_from = source, values_from = intersection, names_prefix = "genes_", values_fill =  list(intersection = NA_character_))

# Combine both
go_combined <- left_join(intsize_wide, genes_wide, by = c("term_id", "term_name"))

go_combined

###########################
# http://revigo.irb.hr/
# default + tiny
revigo_cleaned <- clean_revigo("Revigo_BP_Table_pe-genes.tsv")

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
  filter(merged_term_size < 600 & merged_intersection_size_pe_genes >10)


# Plot with ggpattern
library(ggplot2)

ggplot(go_wide_final, aes(x = reorder(merged_name, percentage_pe_genes), y = percentage_pe_genes)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", percentage_pe_genes)), 
            hjust = -0.1, size = 4) +
  coord_flip() +
  labs(
    title = "GO Term Enrichment",
    x = "GO Term",
    y = "Percentage of PE Genes"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center title
  ylim(0, max(go_wide_final$percentage_pe_genes) * 1.1)  # Give space for text above bars

ggsave("07.02_GO_terms_PE-genes.pdf", width = 40, height = 20, units = "cm")
