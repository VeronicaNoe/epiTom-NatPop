suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(eulerr)
  library(dplyr)
  library(data.table)
})
# setwd("/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/data")
# outPlot<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/plots/"
# outDir<-"/home/IPS2/vibanez/Desktop/Q-lab/IBANEZ_etal_2024/Fig1/results/"
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig1/ab_global-methylation")
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/plots/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/results/"

general_indir<-"/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection"
input<-list.files(path=general_indir,pattern = "DMR.methylation", full.names = FALSE)

out<-c()
for(i in input){
  toAnno<-gsub(".methylation", "",i)
  toAnno<-unlist(strsplit(toAnno, ".", fixed = TRUE))
  DMR<-toAnno[2]
  #chr<-unlist(strsplit(toAnno[1], "_", fixed = TRUE))[1]
  anno<-unlist(strsplit(toAnno[1], "_", fixed = TRUE))[1]
  df<- data.table::fread(paste0(general_indir, "/",i), 
                         sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  if(ncol(df)==0){
    nDMR<-0
    tmp<-cbind.data.frame(DMR, anno, nDMR)
    out<-rbind.data.frame(out, tmp)
  }else{
    nDMR<-nrow(df)
    tmp<-cbind.data.frame(DMR, anno, nDMR)
    out<-rbind.data.frame(out, tmp)
  }
}
head(out)

# Clean chr_size
colnames(chr_size) <- c("chr", "size")
# Total genome size
total_genome_size <- sum(chr_size$size)

# 1. Add total_bp per annotation (each DMR is 100 bp)
out <- out %>%
  mutate(total_bp = nDMR * 100)

# 2. Get total number of DMRs per DMR type (C-DMR or CG-DMR)
total_nDMR_per_DMRtype <- out %>%
  group_by(DMR) %>%
  summarise(total_nDMR = sum(nDMR), 
            total_bp_dmr = sum(total_bp), .groups = "drop")

# 3. Add the genome total bp covered by DMRs (C- + CG- DMRs)
total_bp_all_dmrs <- sum(total_nDMR_per_DMRtype$total_bp_dmr)

# 4. Merge back
out <- out %>%
  left_join(total_nDMR_per_DMRtype, by = "DMR")

# 5. Add columns
out <- out %>%
  mutate(
    genome_coverage_fraction = total_bp_dmr / total_genome_size,
    genome_coverage_percent = genome_coverage_fraction * 100,
    dmr_fraction_within_type = nDMR / total_nDMR,  # proportion within DMR type
    dmr_percent_within_type = dmr_fraction_within_type * 100
  )
#
# Sum coverage for each DMR type
generalCove <- out %>%
  group_by(DMR) %>%
  summarize(covDMR = unique(genome_coverage_percent)) %>%
  ungroup()

# Add Total line manually
total_coverage <- sum(generalCove$covDMR)

generalCove <- generalCove %>%
  add_row(DMR = "Total", covDMR = total_coverage)

# Reorder manually if you want
generalCove$DMR <- factor(generalCove$DMR, levels = c("Total", "C-DMR", "CG-DMR"))

# View
generalCove

ggplot(generalCove, aes(x=DMR, y=covDMR, fill=DMR)) +
  scale_x_discrete(limits = c("Total","C-DMR","CG-DMR"))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("C-DMR"="#820a86", "CG-DMR"="#ffca7b", "Total"="#999999"))+
  ylab("DMR genome coverage (%)") +
  theme_minimal()
ggsave(paste0(outPlot,"05.01.DMR_genomeCoverage.pdf"), width = 60, height = 40, units = "cm")


ggplot(out, aes(x=DMR, y=dmr_percent_within_type, fill=anno)) +
  scale_x_discrete(limits = c("C-DMR","CG-DMR"))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#ccbce2","#9d7ec7", "#d5d3d4","#252123","#b68090", "#f99299"))+
  theme_minimal()
ggsave(paste0(outPlot,"05.02.DMR_featureCoverage.pdf"), width = 60, height = 40, units = "cm")


ggplot(out, aes(x = "", y = dmr_percent_within_type, fill = anno)) +
  geom_bar(stat = "identity", width = 1) +  
  geom_text(aes(label = round(dmr_percent_within_type )), position = position_stack(vjust = 0.5)) +  # Add labels inside slices
  coord_polar("y", start = 0) +  # Make the chart polar
  facet_wrap(~DMR) +  # Facet by DMR
  scale_fill_manual(values = c("#ccbce2", "#9d7ec7", "#d5d3d4", "#252123", "#b68090", "#f99299")) +
  theme_void()  # Remove unnecessary elements
ggsave(paste0(outPlot,"05.03.DMR_featureCoverage_pie.pdf"), width = 60, height = 40, units = "cm")

## venn for poly
vennColors<-c("CG-DMR" = "#820a86","C-DMR"= "#ffca7b", "CG-DMR&C-DMR"="#c77480")
polyepialleles <- euler(c("CG-DMR" = 293047, "C-DMR" = 1794110, "CG-DMR&C-DMR" = 76773), 
                        shape = "ellipse", quantity = quantities)
pdf(paste0(outPlot, "05.04.DMR_overlap_venn.pdf"))
plot(polyepialleles,
     fills= c("#ffca7b","#820a86","#c77480"),
     quantities = list(type = c("percent", "counts"), font = 3))
dev.off()


#=================== GENE LANDSCAPE
dataDir<-"landscape_data/"
input<-list.files(path=dataDir,pattern = "DMR.methylation", full.names = FALSE)

DMRcounts<-data.table::fread("DMR-count.tmp", sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
head(DMRcounts)
rownames(DMRcounts)<-DMRcounts[,1]
out<-c()
for(i in input){
  df<- data.table::fread(paste0(dataDir,i), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
  if(ncol(df)==0){
    nDMR<-0
    tmp<-cbind.data.frame(chr, DMR, anno, nDMR)
    out<-rbind.data.frame(out, tmp)
    next
  }else{
    # save feature names
    toAnno<-gsub(".methylation", "",i)
    print(toAnno)
    toAnno<-unlist(strsplit(toAnno, ".", fixed = TRUE))
    DMR<-toAnno[2]
    chr<-unlist(strsplit(toAnno[1], "_", fixed = TRUE))[1]
    anno<-unlist(strsplit(toAnno[1], "_", fixed = TRUE))[2]
    if (anno == "exon" | anno == "intron"){
      pos<-anno
      anno<-"gene"
    } else {
      pos<-unlist(strsplit(toAnno[1], "_", fixed = TRUE))[3]
      pos<-unlist(strsplit(pos, ".", fixed = TRUE))[1]
    }
    if(anno=="TE-wo-gene" & is.na(pos)){
      pos<-anno
      anno<-anno
    }
    nDMR<-nrow(df)
    totalDMR<-DMRcounts[paste0(chr,"_",DMR),"V2"]
    tmp<-cbind.data.frame(chr, DMR, anno, pos,nDMR, totalDMR)
    out<-rbind.data.frame(out, tmp)
  }
}
out[1:60,]
out$percDMR<-(out$nDMR/out$totalDMR)*100
head(out)
summGene<-plyr::ddply(out, c("DMR", "anno","pos"), summarise,
                   mean = mean(percDMR, na.rm=TRUE),
                   sem = sd(percDMR,na.rm=TRUE)/sqrt(length(percDMR)))

up<-c("arriba-2500","arriba-1500","arriba-1000","arriba-0500")
down<-c("down-0500","down-1000","down-1500","down-2500")

gene<-subset(summGene, anno=="gene")
ggplot(gene, aes(x=pos, y=mean, fill=DMR)) + 
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(1))+
  scale_x_discrete(limits = c(up,"exon",'intron',down))+
  scale_fill_manual(values=c("#820a86", "#ffca7b"))+
  facet_grid(.~anno)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("DMR percentage")
ggsave(paste0(outPlot,"05.05_gene-DMR_landscape.pdf"), width = 60, height = 40, units = "cm")

te<-subset(summGene, anno=="TE-wo-gene")
ggplot(te, aes(x=pos, y=mean, fill=DMR)) + 
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(1))+
  scale_x_discrete(limits = c(up,"TE-wo-gene",down))+
  scale_fill_manual(values=c("#820a86", "#ffca7b"))+
  facet_grid(.~anno)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("DMR percentage")
ggsave(paste0(outPlot,"05.06_TE-DMR_landscape.pdf"), width = 60, height = 40, units = "cm")
