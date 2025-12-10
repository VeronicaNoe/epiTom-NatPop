suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
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

ko_targets_file<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/results/07.1_DMRs-over-KO.tsv"
ko_targets<-fread(ko_targets_file, sep = '\t', na.strings = "NA", header = T)

#ADD GENES

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
axis_df <- all_dt[, .(center = (min(BPcum) + max(BPcum)) / 2), by = CHR]
all_dt_filtered <- all_dt[logP >= 2]

ggplot(all_dt_filtered, aes(x = BPcum, y = logP)) +
  geom_point(aes(shape = molmarker, fill = molmarker, color = factor(CHR)),
             size = 2, alpha = 0.8) +
  scale_shape_manual(values = c(SNP = 19, DMR = 25)) +
  scale_fill_manual(values = c(SNP = "#999999", DMR = "#c06a80")) +
  scale_color_manual(values = rep(c("black", "grey60"), length(unique(all_dt$CHR)))) +
  geom_hline(yintercept = -log10(snp_thresh), linetype = "dashed", color = "#999999") +
  geom_hline(yintercept = -log10(dmr_thresh), linetype = "dashed", color = "#c06a80") +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  theme_bw() +
  labs(title = paste0(meta," Manhattan Plot"), x = "Chromosome", y = "-log10(P)") +
  theme(
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(out_plot, meta,"_combined_manhattan.pdf"), width = 12, height = 5, dpi = 300)

###
# Create a single data.table with marker type
qq_data <- rbind(
  data.table(P = snp_dt$P, molmarker = "SNP"),
  data.table(P = dmr_dt$P, molmarker = "DMR")
)
# Remove NAs and zeroes
qq_data <- qq_data[!is.na(P) & P > 0]
# Compute expected and observed -log10(P)
qq_data[, expected := -log10(ppoints(.N)), by = molmarker]
qq_data[, observed := -log10(sort(P)), by = molmarker]
# Function to compute genomic inflation factor
calc_lambda <- function(pvals) {
  chisq <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  return(round(lambda, 3))
}

lambda_snp <- calc_lambda(snp_dt$P)
lambda_dmr <- calc_lambda(dmr_dt$P)
qq_title <- paste0("QQ Plot - ",meta," (SNP vs DMR)")

# Plot
ggplot(qq_data, aes(x = expected, y = observed, color = molmarker, shape = molmarker, fill = molmarker)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  scale_color_manual(values = c(DMR = "#c06a80", SNP = "#999999")) +
  scale_fill_manual(values = c(DMR = "#c06a80", SNP = "#999999")) +
  scale_shape_manual(values = c(SNP = 16, DMR = 25)) +
  theme_bw() +
  labs(title = qq_title,
       x = "Expected -log10(P)", y = "Observed -log10(P)") +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+ 
  annotate("text", x = 1, y = max(qq_data$observed), hjust = 0,
           label = sprintf("λ_SNP = %.3f\nλ_DMR = %.3f", lambda_snp, lambda_dmr),
           size = 4, color = "black")

ggsave(paste0(out_plot, meta,"_combined_qqplot.pdf"), width = 5.5, height = 5.5, dpi = 300)

}
