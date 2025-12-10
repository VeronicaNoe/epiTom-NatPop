suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
})
# Input and output directories
inDir<-"/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection/bb_pedigree/"
pedResult<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/af_pedigree/"
pedPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/plots/"

dmrTable<-fread(paste0(pedResult,"05.0_DMRs-pedigree.tsv"), sep = '\t', header = T, 
                data.table = T, fill = TRUE, na.strings = c(""), nThread = 20)
dmrTable$pos<-paste(dmrTable$V1,dmrTable$V2, sep="_")
input<-list.files(path=inDir, pattern = '.methylation', full.names = FALSE)
input<-grep('all', input, value=T)
df<-c()
for (i in input){
  dmr<-unlist(strsplit(i, '_', fixed = T))[1]
  comparison<-unlist(strsplit(i, '_', fixed = T))[2]
  annotation<-unlist(strsplit(i, '_', fixed = T))[3]
  annotation<-unlist(strsplit(annotation, '.', fixed = T))[1]
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
  tmp$pos<-paste(tmp$V1,tmp$V2, sep="_")
  tmp$DMR<-dmr
  tmp<-tmp %>%
    filter(`DMR`==dmr & pos %in% dmrTable$pos)
  if (is.null(tmp) ) {
    next
  }
  tmp$annotation<-annotation
  dmrTmp<-dmrTable %>%
    filter(pos %in% tmp$pos)
  
  df<-rbind.data.frame(df,tmp)
}

data.table::fwrite(df, file = paste0(pedResult,"05.3_DMRs-pedigree_annotated.tsv"), quote = FALSE, 
                   row.names = FALSE, col.names = T, sep = "\t")
toRemove<-c("3UTR","5UTR","exon","intron", "intron-TE", "exon-TE")
df<-df%>%
  filter(!annotation %in% toRemove)

summary <- df %>%
  group_by(DMR) %>%
  mutate(total_per_group = n()) %>% # Total rows per comparison and DMR
  group_by(DMR, annotation) %>%
  summarise(
    total_per_group = first(total_per_group), # Carry forward total rows per group
    count_per_annotation = n(),              # Count of rows per annotation
    percentage = (count_per_annotation / total_per_group) * 100, # Percentage calculation
    .groups = "drop"
  )

#### get info per comparison
input<-list.files(path=inDir, pattern = '.methylation', full.names = FALSE)
input<-grep('all', input, value=T, invert = T)
perCom<-c()
for (i in input){
  dmr<-unlist(strsplit(i, '_', fixed = T))[1]
  comparison<-unlist(strsplit(i, '_', fixed = T))[2]
  ann<-unlist(strsplit(i, '_', fixed = T))[3]
  ann<-unlist(strsplit(ann, '.', fixed = T))[1]
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
  tmp$pos<-paste(tmp$V1,tmp$V2, sep="_")
  tmp$DMR<-dmr
  tmpAnno<-df %>%
    filter(annotation == ann & DMR==dmr)
  tmp<-tmp %>%
    filter(pos %in% tmpAnno$pos)
  if (is.null(tmp) ) {
    next
  }
  tmp$annotation<-ann
  counts<-nrow(tmp)
  totalAnnotation<-nrow(tmpAnno)
  tmpOut<-cbind(dmr, ann, comparison,totalAnnotation, counts)
  perCom<-rbind.data.frame(perCom,tmpOut)
}


sumPerCom<-perCom %>%
  filter(!ann %in% toRemove) %>%
  group_by(dmr, comparison, ann) %>%
  summarise(
    total_counts = sum(as.numeric(counts), na.rm = TRUE), # Sum counts per group
    .groups = "drop"                         # Remove grouping after summarization
  ) %>%
  group_by(dmr, comparison) %>%
  mutate(
    total_per_group = sum(total_counts),  # Total counts per dmr and comparison
    proportion = total_counts / total_per_group * 100 # Proportion per annotation
  ) %>%
  ungroup() # R

pie_chart <- ggplot(sumPerCom, aes(x = "", y = proportion, fill = ann)) +
  geom_bar(stat = "identity", width = 1) +  # Create bar plot
  geom_text(
    aes(label = paste0(round(proportion, 1), "%")),  # Label with rounded proportions
    position = position_stack(vjust = 0.5),  # Center the text in the slices
    color = "white",  # White text
    size = 3  # Adjust text size
  ) +
  coord_polar(theta = "y") +  # Convert to pie chart
  facet_wrap(~ dmr + comparison, scales = "free") +  # Facet by dmr and comparison
  labs(
    title = "Proportion of Annotations per DMR and Comparison",
    x = NULL,
    y = NULL,
    fill = "Annotation"
  ) +
  scale_fill_manual(values = c("#ccbce2", "#9d7ec7", "#d5d3d4", "#252123", "#b68090")) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    panel.grid = element_blank()  # Remove grid
  )

# Save the updated pie chart
ggsave(paste0(pedPlot, "05.4_pie_with_labels.pdf"), plot = pie_chart, width = 8, height = 6)
