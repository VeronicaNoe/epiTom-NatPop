suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(stringr)
  library(ggrepel)
})

# === FILES ===
in_dir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results"
out_plot<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/plots/sigTraits_combined_manhatan/"
ld_file<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/ab_get-LD/bb_general/C-DMR_LD.ld.gz"
ld<-fread(ld_file, sep = ' ', na.strings = "NA", header = T)
ld[, CHR_A := gsub("SL2.50", "", CHR_A)]
ld[, CHR_B := gsub("SL2.50", "", CHR_B)]
ld[, SNP_A := gsub("SL2.50", "", SNP_A)]
ld[, SNP_B := gsub("SL2.50", "", SNP_B)]


#ADD GENES
gene_dir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/results/"
input<-list.files(path=gene_dir, pattern = "closest", full.names = FALSE)
input<-grep('04', input, value=T)
out<-c()
for(i in input){
  file_parts <- unlist(strsplit(i, "\\_"))  # Split on period
  ann <-gsub('-closest-','_',file_parts[2])
  dt<-data.table::fread(paste0(gene_dir,i),sep = '\t', data.table = T, header = F, 
                        fill = TRUE, na.string=c("NA"), nThread = 20)
  dt<-dt %>% 
    filter(V5!="SNP", V7=="aBN", V8=="sig")
  if(ann=="QTL"){
    annotation<-"QTL"
    toKeep<-c(4,9:12,14)
  }else{
    annotation<-"QTLoverTE"
    toKeep<-c(4,9,16:18,20)
  }
  dt<-dt[,..toKeep]
  colnames(dt)<-c('metabolite','QTL','chr','start','end', 
                  'geneName')
  dt$annotation<-annotation
  out<-rbind.data.frame(out,dt)
}
out$geneName <- gsub("Name=","",out$geneName)
out <- out %>%
  mutate(geneName = str_replace(geneName, "\\.\\d+$", ""))

nearby_genes<-data.table::fread(paste0(gene_dir,"04.3_nearby-genes.tsv"),
                                sep = '\t', data.table = T, header = F, 
                                fill = TRUE, na.string=c("NA"), nThread = 20)

nearby_genes <- nearby_genes %>%
  filter(V5!="SNP", V7=="aBN", V8=="sig") %>%
  select(V4,V9,V10,V11,V12,V14) %>%
  mutate(geneName = gsub("Name=","",V14),
         geneName = str_replace(geneName, "\\.\\d+$", ""),
         annotation = "nearby-closest_to_top_mqtl")%>%
  select(-V14)
colnames(nearby_genes)<-c('metabolite','QTL','chr','start','end', 
                'geneName', 'annotation')
# Join to get nearby genes matched by geneName
closest_genes <- rbind(out, nearby_genes)

#=====
trait_info<-fread("/mnt/disk2/vibanez/10_data-analysis/Fig5/results/01.0_general-information.tsv", 
                  sep = '\t', na.strings = "NA", header = T)
# traits<-trait_info%>%
#   filter(mm=="DMR", resultType=="sig",kinship=="DMR-SNP", index=="aBN" ) %>%
#   pull(unique(metabolite))

#traits<-c("Tm283","Tm220","Tm242","Tm205")
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
    snps <- chr_dt[molmarker == "SNP"]
    snps[, R2 := NA_real_]  # or R2 := 0 if you prefer zeros
    
    # Get top DMR (assume only one has is_sig == TRUE)
    top_dmr <- chr_dt[is_sig == TRUE & molmarker == "DMR"]
    top_id <-sub(":[A-Z]+-DMR$", "", top_dmr$V1)
    ld_with_top <- ld[SNP_A == top_id | SNP_B == top_id]
    
    # STEP 3: For convenience, always put the other DMR in a column `partner`
    ld_with_top[, partner := ifelse(SNP_A == top_id, SNP_B, SNP_A)]
    # STEP 3: For convenience, always put the other DMR in a column `partner`
    ld_with_top[, partner := ifelse(SNP_A == top_id, SNP_B, SNP_A)]
    ld_with_top[, partner_dmr := paste0(partner, ":C-DMR")]
    
    # STEP 5: Join DMRs with LD R2 to the top DMR
    dmrs_chr <- merge(dmrs,
                      ld_with_top[, .(partner_dmr, R2)],
                      by.x = "SNP",by.y = "partner_dmr",all.x = TRUE)
    
    dmrs_chr[is.na(R2), R2 := 0]
    
    chr_dt <- rbind(dmrs_chr, snps)
    
    #
    # Filter genes to the current region (chr and window)
    genes_in_region <- closest_genes[
      metabolite == meta & chr == chr_target &
        start >= bp_min & end <= bp_max]
    
    # Shift y-position below the Manhattan plot
    gene_y <- -0.5
    
    ggplot(chr_dt, aes(x = BP, y = logP)) +
      # SNPs - grey points
      geom_point(
        data = chr_dt[molmarker == "SNP" & is_sig == FALSE],
        aes(shape = molmarker),
        color = "#999999", fill = "#999999",
        size = 2, alpha = 0.4
      ) +
      geom_point(
        data = chr_dt[molmarker == "SNP" & is_sig == TRUE],
        aes(shape = molmarker),
        color = "#999999", fill = "#999999",
        size = 3.5, stroke = 1.1
      ) +
      
      # DMRs - colored by R²
      geom_point(
        data = chr_dt[molmarker == "DMR" & is_sig == FALSE],
        aes(shape = molmarker, color = R2, fill = R2),
        size = 2, alpha = 0.4
      ) +
      geom_point(
        data = chr_dt[molmarker == "DMR" & is_sig == TRUE],
        aes(shape = molmarker, color = R2, fill = R2),
        size = 3.5, stroke = 1.1
      ) +
      
      # Gene bars
      geom_segment(
        data = genes_in_region,
        aes(x = start, xend = end, y = gene_y, yend = gene_y),
        inherit.aes = FALSE,
        color = "darkblue", size = 1
      ) +
      
      # Gene names with ggrepel
      geom_text_repel(
        data = genes_in_region,
        aes(x = (start + end) / 2, y = -log10(dmr_thresh) + 0.2, label = geneName),
        inherit.aes = FALSE,
        direction = "y",
        angle = 45,
        size = 2.6,
        segment.color = "grey70",
        min.segment.length = 0,
        box.padding = 0.15,
        max.overlaps = Inf
      ) +
    
      # Highlight top DMR
      geom_vline(
        xintercept = top_dmr$BP,
        linetype = "dashed",
        color = "red",
        size = 0.6
      ) +
      
      # Color scale for LD
      scale_color_viridis_c(option = "A", direction = 1, name = expression(LD~(R^2)), limits = c(0, 0.03)) +
      scale_fill_viridis_c(option = "A", direction = 1, name = expression(LD~(R^2)), limits = c(0, 0.03)) +
      
      # Shapes
      scale_shape_manual(values = c(SNP = 16, DMR = 25)) +
      
      # Significance thresholds
      geom_hline(yintercept = -log10(snp_thresh), linetype = "dashed", color = "#999999") +
      geom_hline(yintercept = -log10(dmr_thresh), linetype = "dashed", color = "#c06a80") +
      
      labs(
        title = paste0(meta, " Manhattan Plot (Chr ", chr_target, ": Mb)"),
        x = "Position (bp)", y = expression(-log[10](P))
      ) +
      
      theme_bw() +
      theme(
        legend.title = element_text(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    ggsave(paste0(out_plot, meta, "_chr", chr_target, "_combined_manhattan_LD",qtl,".pdf"),
           width = 10, height = 4, dpi = 300)
    
    }

  }
