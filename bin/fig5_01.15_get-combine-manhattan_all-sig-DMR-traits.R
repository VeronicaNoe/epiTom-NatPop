suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(tidyr)
})

# === FILES ===
in_dir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results"
out_plot<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/plots/sigTraits_combined_manhatan/"

ko_targets_file<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/results/07.1_DMRs-over-KO.tsv"
ko_targets<-fread(ko_targets_file, sep = '\t', na.strings = "NA", header = T)

#ADD GENES

#=====
trait_info<-fread("/mnt/disk2/vibanez/10_data-analysis/Fig5/results/01.0_general-information.tsv", 
                  sep = '\t', na.strings = "NA", header = T)
# traits<-trait_info%>%
#   filter(mm=="DMR", resultType=="sig",kinship=="DMR-SNP", index=="aBN" ) %>%
#   pull(unique(metabolite))
traits<-c("Tm337")
mqtl_info<-fread("/mnt/disk2/vibanez/10_data-analysis/Fig5/results/01.1_QTLs_with_coordinates.tsv", 
                 sep = '\t', na.strings = "NA", header = T)
mqtl_info <-mqtl_info %>%
  filter(mm=="DMR", resultType=="sig",kinship=="DMR-SNP", index=="aBN", metabolite %in% traits)
# === Function to Read and Format ===
read_and_format <- function(file, moltype) {
  dt <- fread(file, sep = '\t', na.strings = "NA", header = FALSE)
  dt[, CHR := as.numeric(sub("ch", "", sapply(strsplit(V1, ":"), `[[`, 1)))]
  dt[, BP := as.numeric(sapply(strsplit(V1, ":"), `[[`, 2))]
  dt[, SNP := V1]
  dt[, P := V4]
  dt[, molmarker := moltype]
  dt[, logP := -log10(P)]
  #dt[is.na(logP), logP := 0]
  dt
}
# === COLORS (qqman style) ===
chr_colors <- rep(c("black", "grey"), 25)
# === PARAMETERS ===
snp_thresh <- 4.1e-08
dmr_thresh <- 1.335748e-07

## 
for(meta in traits){
  snpResult<-trait_info%>%
    filter(metabolite==meta, kinship=="DMR-SNP", index=="aBN", mm=="SNP") %>%
    pull(resultType)
  snp_file <- paste0(in_dir,"/",snpResult,"/",meta,".leaf-metabolites_SNP_DMR-SNP_aBN.ps.gz")
  dmrResult<-trait_info%>%
    filter(metabolite==meta, kinship=="DMR-SNP", index=="aBN", mm=="DMR") %>%
    pull(resultType)
  dmr_file <- paste0(in_dir,"/",dmrResult,"/",meta,".leaf-metabolites_DMR_DMR-SNP_aBN.ps.gz")
  # === Read Data ===
  snp_dt <- read_and_format(snp_file, "SNP")
  dmr_dt <- read_and_format(dmr_file, "DMR")
  all_dt <- rbind(snp_dt, dmr_dt)
  # ====
  # get chromosome sizes
  chr_sizes <- all_dt[, .(chr_len = max(BP)), by = CHR][order(CHR)]
  chr_sizes[, chr_start := cumsum(shift(chr_len, fill = 0))]
  # Merge cumulative position
  all_dt <- merge(all_dt, chr_sizes[, .(CHR, chr_start)], by = "CHR")
  all_dt[, BPcum := BP + chr_start]
  # === Filter to Chromosome in Region ===
  sub_mqtl<-mqtl_info %>%
    filter(metabolite==meta) 
  
  # === PARAMETERS ===
  for(i in 1:nrow(sub_mqtl)){
    qtl<-gsub(":","_",as.character(sub_mqtl[i,"qtl"]))
    chr_target <-as.numeric(sub_mqtl[i,"chr"])
    bp_min <- as.numeric(sub_mqtl[i,"start"]-500000)
    bp_max <- as.numeric(sub_mqtl[i,"end"]+500000)
    chr_dt <- all_dt[CHR == chr_target & BP >= bp_min & BP <= bp_max]
    chr_dt[, is_sig := (molmarker == "SNP" & P <= snp_thresh) |
             (molmarker == "DMR" & P <= dmr_thresh)]
    # Get all DMRs
    dmrs <- chr_dt[molmarker == "DMR"]
    # STEP 1: Prepare matching DMR IDs from ko_targets to match V1
    # Convert "ch02:25695001:C-DMR" → "ch02_C-DMR_25695001"
    dmrs[, dmr_id := paste0(gsub(":", "", sub(":.*", "", V1)), "_C-DMR_", BP)]
    dmrs[, DMR_order := .I]  # keep original row order
    # STEP 2: Filter ko_targets to only the DMRs found in your `dmrs` table
    # Set bin size
    bin_size <- 20000
    ko_sub <- ko_targets[DMR %in% dmrs$dmr_id]
    # Extract numeric position from DMR ID
    ko_sub[, chr := gsub("SL2.50", "", chr)]  # optional cleanup
    ko_sub[, bin := floor(start / bin_size) * bin_size]
    # Split KO string into one row per KO
    ko_long <- ko_sub[, .(ko = unlist(strsplit(ko, ","))), by = .(DMR, chr, start, bin)]
    heatmap_df <- ko_long[, .N, by = .(ko, chr, bin)]
    
    # ...
    ggplot(heatmap_df, aes(x = bin, y = ko, fill = N)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "C", direction = -1, name = "# DMRs") +
      labs(title = "DMR_KO Density per Genomic Region",
        x = "Genomic Position (binned, bp)",
        y = "Knockout") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    ggsave(paste0(out_plot, meta, "_chr", chr_target, "_ko-targets",qtl,".pdf"), width = 10, height = 4, dpi = 300)
  }
}
