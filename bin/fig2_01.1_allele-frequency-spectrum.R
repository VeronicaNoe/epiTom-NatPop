suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(tidyr)
})
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig2/aa_get-allele-spectrum")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/plots/"

input<-list.files(pattern = ".edited-count", full.names = FALSE)
sfs_cum <- list()
sfs_histo <- list()
for (s in input) {
    print(s)
    # Extract annotation and marker information
    group <- unlist(strsplit(s, "\\_"))[1]
    anno <- unlist(strsplit(s, "\\_"))[2]
    marker <- unlist(strsplit(s, "\\."))[1]
    marker <- unlist(strsplit(marker, "\\_"))[3]
    # Read the data
    dt <- fread(s, sep = '\t', data.table = TRUE, header = FALSE, nThread = 20, fill = TRUE)
    setnames(dt, c('chr', 'start', 'nAllele', 'nMinorAllele'))
    setorder(dt, chr, start)
    # Filter and mutate using data.table
    dt <- dt[chr != "SL2.50ch00"]
    dt[, MAF := nMinorAllele / nAllele]
    dt <- dt[!is.nan(MAF)]
    # Create histogram and cumulative data
    histo <- hist(dt$MAF, breaks = seq(from = 0, to = 0.5, by = 0.01), plot = FALSE)
    alleleCum <- data.table(
      nAlleles = histo$counts,
      MAF = histo$breaks[-1],  # Use -1 to get the right edge of each bin
      cumulatedAlleles = cumsum(histo$counts) / sum(histo$counts),
      group = group,
      marker = marker,
      annotation = anno)
    # Create a data frame with the relevant data
    hist_table <- data.frame(
      breaks_start = histo$mids,
      breaks_end = histo$breaks[-1],
      counts = histo$counts,
      density = histo$density/100,
      group = group,
      marker = marker,
      annotation = anno
    )
    hist_table <- hist_table %>%
      pivot_longer(cols = c("breaks_start", "breaks_end"), 
                   names_to = "type_breaks", values_to = "MAF")
    # Append the result to the list
    sfs_cum[[length(sfs_cum) + 1]] <- alleleCum
    sfs_histo[[length(sfs_histo) + 1]] <- hist_table
}
# Combine all results into a single data.table
sfsCum <- rbindlist(sfs_cum)
sfsHisto <- rbindlist(sfs_histo)
data.table::fwrite(sfsCum, "01.01.SFS_cumulated.tsv", quote=F, row.names=F,sep="\t",nThread=12)
data.table::fwrite(sfsHisto, "01.01.SFS_histogram-values.tsv", quote=F, row.names=F,sep="\t",nThread=12)


generalCum<-sfsCum %>%
  filter(group=="allAcc")

generalHist<-sfsHisto %>%
  filter(group=="allAcc")

for(a in unique(generalCum$annotation)){
  tmp<-generalCum %>%
    filter(annotation==a)
  ggplot(tmp, aes(x=MAF, y=cumulatedAlleles, color=marker)) +
    geom_line(aes(), stat = "identity")+
    geom_line(aes(x=0.05), color="red", linetype="dashed")+
    annotate("text", x = 0.030, y = 0.95, angle = 90, label = "MAF 0.05", 
             vjust = 0.5) +
    scale_x_continuous(breaks = seq(0.00,0.5, by=0.05))+
    scale_y_continuous(breaks = seq(0,1, by=0.1))+
    scale_color_manual(values=c("#820a86", "#ffca7b", "#999999"))+
    xlab("Minor Allele Frequency") + ylab("Cumulated number of observed alleles")+
    ggtitle(paste0("MAF for ", a, " annotation")) #+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))  # Rotate x-axis labels horizontally
  ggsave(paste(outPlot,"01.01.all-accessions_",a,"_SFS-cumulated-normalized.pdf"),  width = 20, height = 20, units = "cm")
  
  tmpHisto<-generalHist %>%
    filter(annotation==a)
  ggplot(tmpHisto, aes(x=MAF, y=density, color=marker)) +
    geom_line(aes(), stat = "identity")+
    geom_line(aes(x=0.05), color="red", linetype="dashed")+
    annotate("text", x = 0.030, y = 0.40, angle = 90, label = "MAF 0.05", 
             vjust = 0.5) +
    scale_x_continuous(breaks = seq(0.00,0.5, by=0.05))+
    #scale_y_continuous(breaks = seq(0,0.5, by=0.1))+
    scale_color_manual(values=c("#820a86", "#ffca7b", "#999999"))+
    xlab("Minor Allele Frequency") + ylab("Allele density")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(paste0("MAF for ", a, " annotation")) +
    theme_bw()
  ggsave(paste0(outPlot,"01.02.all-accessions_",a,"_SFS-histogram-normalized.pdf"),  width = 20, height = 20, units = "cm")
}

## by groups
groupsCum<-sfsCum %>%
  filter(group!="allAcc")

groupsHist<-sfsHisto %>%
  filter(group!="allAcc")

ggplot(groupsCum, aes(x=MAF, y=cumulatedAlleles, color=marker)) +
    geom_line(aes(), stat = "identity")+
    geom_line(aes(x=0.05), color="red", linetype="dashed")+
    annotate("text", x = 0.030, y = 0.95, angle = 90, label = "MAF 0.05", 
             vjust = 0.5) +
    facet_grid(group~., scales = "free")+
    scale_x_continuous(breaks = seq(0.00,0.5, by=0.05))+
    scale_y_continuous(breaks = seq(0,1, by=0.1))+
    scale_color_manual(values=c("#820a86", "#ffca7b", "#999999"))+
    xlab("Minor Allele Frequency") + ylab("Cumulated number of observed alleles")+
    ggtitle("MAF for general annotation") #+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))  # Rotate x-axis labels horizontally
  ggsave(paste(outPlot,"01.03.byGroups_SFS-cumulated-normalized.pdf"),  width = 20, height = 20, units = "cm")
  
  ggplot(groupsHist, aes(x=MAF, y=density, color=marker)) +
    geom_line(aes(), stat = "identity")+
    geom_line(aes(x=0.05), color="red", linetype="dashed")+
    facet_grid(group~., scales = "free")+
    annotate("text", x = 0.030, y = 0.40, angle = 90, label = "MAF 0.05", 
             vjust = 0.5) +
    scale_x_continuous(breaks = seq(0.00,0.5, by=0.05))+
    #scale_y_continuous(breaks = seq(0,0.5, by=0.1))+
    scale_color_manual(values=c("#820a86", "#ffca7b", "#999999"))+
    xlab("Minor Allele Frequency") + ylab("Allele density")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("MAF for general annotation") +
    theme_bw()
  ggsave(paste0(outPlot,"01.04.byGroups_SFS-histogram-normalized.pdf"),  width = 20, height = 20, units = "cm")

# dim(sfs)
# subDF<-subset(sfs, MAF>=0.05 & marker!="SNPs")
# dim(subDF)
# data.table::fwrite(subDF, paste0(outDir,"DMR_filteredByMAF.tsv"),quote=F,row.names=F,sep="\t",nThread=12)
# subDF<-subset(sfs, marker!="SNPs")
# dim(subDF)
# data.table::fwrite(subDF, paste0(outDir,"DMR_MAFs.tsv"),quote=F,row.names=F,sep="\t",nThread=12)
# 
# #make a list of DMR with MAF>0.05
# out<-subset(sfs, MAF>=0.05)
# toList<-as.data.frame(paste0(out$V1,'_',out$marker,'_',out$V2,'.phenotype'))
# data.table::fwrite(toList, paste0(outDir,"toTargetMAF005.tsv"),quote=F,row.names=F,sep="\t",nThread=12)
# # then, tho this in the terminal
# #cat toTargetMAF005.tsv |awk '{OFS="\t"}NR>1{print}' | sed 's/SL2.50//g' | sed 's/.edited-count//g' > toTarget-MAF005
# 
