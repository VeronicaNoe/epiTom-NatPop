#!/users/bioinfo/vibanez/miniconda3/envs/methylKit/bin/Rscript
require("data.table")
library(ggplot2)
library(plyr)
setwd("/mnt/disk2/vibanez/02_methylkit/ak_epiPI/aa_data")
outDir<-"/mnt/disk2/vibanez/02_methylkit/ak_epiPI/ab_results/"
outPlot<-"/mnt/disk2/vibanez/02_methylkit/ak_epiPI/ac_plots/"
# by features

infoCols<-c("chr", "start", "end")
{
#input<-list.files(pattern = ".bed", full.names = FALSE)
anno<-c("gene","gene-TE","intergenic","promotor", "TE-wo-gene", "general")
diversity<-c()
for(a in anno){
  pi<- data.table::fread(paste0('SNPs_',a,"_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  pi$V8<-NULL
  colnames(pi)<-c(infoCols,'1-WILD_SNPs','2-PIM_SNPs','3-SLC_SNPs','4-SLL_SNPs')
  cDMR<- data.table::fread(paste0('C-DMR_',a,"_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  cDMR$V8<-NULL
  colnames(cDMR)<-c(infoCols,'1-WILD_C-DMR','2-PIM_C-DMR','3-SLC_C-DMR','4-SLL_C-DMR')
  cgDMR<- data.table::fread(paste0('CG-DMR_',a,"_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
  cgDMR$V8<-NULL
  colnames(cgDMR)<-c(infoCols,'1-WILD_CG-DMR','2-PIM_CG-DMR','3-SLC_CG-DMR','4-SLL_CG-DMR')
  
  merged_table <- dplyr::left_join(pi, cDMR, by = infoCols)
  merged_table <- dplyr::left_join(merged_table, cgDMR, by = infoCols)
  head(merged_table)
  merged_table$annotation<-a
  diversity<-rbind.data.frame(diversity, merged_table)
}
head(diversity)
diversity <- reshape2::melt(diversity, id.vars = c("chr", "start", "end", "annotation"), variable.name = "samples",
                    value.name = "diversity")
head(diversity)

summ<-plyr::ddply(diversity, c("samples","annotation"), summarise,
                   mean = mean(diversity, na.rm=TRUE)*1000, # 1e3
                   sd= sd(diversity, na.rm=TRUE)*1000,
                   sem = sd(diversity,na.rm=TRUE)/sqrt(length(diversity))*1000)
summ
summ <- tidyr::separate(summ, samples, into = c("group", "marker"), sep = "_", remove = T)


ggplot(summ, aes(x=group, y=mean, fill=group))+
  geom_bar(stat='identity', position = position_dodge())+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(1))+
  scale_fill_manual(values=c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(annotation ~marker, scales = 'free_x')+
  coord_flip()
ggsave(paste0(outDir,"diversity-by-feature_20acc.png"),  width = 20, height = 20, units = "cm")  

diversity <- tidyr::separate(diversity, samples, into = c("group", "marker"), sep = "_", remove = T)
data.table::fwrite(diversity, file=paste0(outDir,"diversity_complete.tsv"), quote=F,row.names=F,col.names = T,sep="\t")

ggplot(diversity, aes(x = group, y = diversity, color=group, fill=group, alpha=0.4)) +
  geom_violin( na.rm = T) +
  geom_boxplot(width=0.5, na.rm = T, alpha=0.1,lwd = 0.3,outlier.shape = NA, 
	color= c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  stat_summary(fun = mean, geom = "point", shape = 16, size = 2.5, color = "darkred", 
	position = position_dodge(width = 0.75)) + 
  labs(x = "", y = "Diversity", title = "") +
	scale_fill_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
	scale_color_manual(values = c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  theme_bw() +
facet_grid(annotation ~marker, scales = 'free_x')+
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))
ggsave(paste0(outPlot,"diversity-by-feature.pdf"),  width = 20, height = 20, units = "cm")  



}
###############
# general
{
pi<- data.table::fread(paste0("SNPs_general_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
head(pi)
pi$V8<-NULL
colnames(pi)<-c(infoCols,'1-WILD_SNPs','2-PIM_SNPs','3-SLC_SNPs','4-SLL_SNPs')

cDMR<- data.table::fread(paste0("C-DMR_general_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
cDMR$V8<-NULL
colnames(cDMR)<-c(infoCols,'1-WILD_C-DMR','2-PIM_C-DMR','3-SLC_C-DMR','4-SLL_C-DMR')
cgDMR<- data.table::fread(paste0("CG-DMR_general_20acc_merged-pi.bed"), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
cgDMR$V8<-NULL
colnames(cgDMR)<-c(infoCols,'1-WILD_CG-DMR','2-PIM_CG-DMR','3-SLC_CG-DMR','4-SLL_CG-DMR')

diversity <- dplyr::left_join(pi, cDMR, by = infoCols)
diversity <- dplyr::left_join(diversity, cgDMR, by = infoCols)
# If you want to replace NA values with 0, you can use the following line
diversity$annotation<-"general"
head(diversity)

piDF <- reshape2::melt(diversity, id.vars = c("chr", "start", "end", "annotation"), variable.name = "samples",
                       value.name = "diversity")
#piDF <- tidyr::separate(piDF, samples, into = c("group", "marker"), sep = "_", remove = T)

split_samples <- strsplit( as.character(piDF$samples), "_", fixed = TRUE)
gc()
piDF$group <- sapply(split_samples, `[`, 1)
gc()
piDF$marker <- sapply(split_samples, `[`, 2)
rm(pi, cDMR, cgDMR, diversity)
gc()
head(piDF)
summ<-plyr::ddply(piDF, c("marker","annotation","group"), summarise,
                  mean = mean(diversity, na.rm=TRUE)*1000, # 1e3
                  sd= sd(diversity, na.rm=TRUE)*1000,
                  sem = sd(diversity,na.rm=TRUE)/sqrt(length(diversity))*1000)

summ
png(paste0(outDir,"/diversity-all-feature_20acc.png"))
ggplot(summ, aes(x=group, y=mean, fill=group))+
  geom_bar(stat='identity', position = position_dodge())+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(1))+
  scale_fill_manual(values=c("#000000ff","#00a681ff","#ed83b5ff", "#00a4deff"))+
  facet_grid(. ~marker, scales = 'free_x')+
  coord_flip()
dev.off()

}
########################## check for correlations 
# check<-subset(summGlobal, DMR =="C-DMR")
# indNumber<-c(17,28,7,30,5,25,17,12,5,9,4,5)
# check$nIndividuals<-indNumber
# plot(x=check$nIndividuals, y=check$mean)
#####################################################
# subGlobal<-subset(summGlobal, GROUP!="3-SLC-wIntro" & 
#                     GROUP!="3-SLC-woIntro" & 
#                     GROUP!="4-SLL-wIntro" & 
#                     GROUP!="4-SLL-woIntro")
# dmr<-unique(subGlobal$DMR)
# for (d in dmr){
#   tmp<-subset(subGlobal, DMR==d)
#   ggplot(tmp, aes(x=GROUP, y=mean, fill=GROUP))+
#     geom_bar(stat='identity', position = position_dodge())+
#     geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
#                   position=position_dodge(1))+
#     scale_fill_manual(values=c( "#209f68","#f6a29a","#d75e5e","#c33131"))+
#     #facet_wrap(~DMR)+
#     coord_flip()+ 
#     theme_minimal()
#   ggsave(paste0(outDir,"1_",d,"_epi_global.png"),  width = 20, height = 20, units = "cm")
# }

################## plot for each annotation
specific<-subset(epiPi, ANNOTATION!="general")
summSpecific<-plyr::ddply(specific, c("GROUP","DMR", "ANNOTATION"), summarise,
                        mean = mean(PI, na.rm=TRUE)*1000, # 1e3
                        sd= sd(PI, na.rm=TRUE)*1000,
                        sem = sd(PI,na.rm=TRUE)/sqrt(length(PI))*1000)

subSpecific<-subset(summSpecific, GROUP!="3-SLC-wIntro" & 
                    GROUP!="3-SLC-woIntro" & 
                    GROUP!="4-SLL-wIntro" & 
                    GROUP!="4-SLL-woIntro")

dmr<-unique(subSpecific$DMR)
ann<-unique(subSpecific$ANNOTATION)

#line plot BUT try it with all dots
ggplot(summSpecific, aes(x=GROUP, y=mean, group=ANNOTATION)) +
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.05,) +
    geom_line(aes(color=ANNOTATION), linewidth=1)+
    geom_point(aes(color=ANNOTATION))+
    scale_color_manual(values = c("#ccbce2","#9d7ec7", "#d5d3d4","#b68090", "#f99299"))+
    facet_grid(. ~DMR, scales = 'free_x')+
    coord_flip()
ggsave(paste0(outDir,"2_feature-specific-epi.png"),  width = 20, height = 20, units = "cm")


# #pie chart
# subSpecific
# group<-unique(subSpecific$GROUP)
# toPie<-c()
# for (g in group){
#   for (d in dmr){
#     subS<-subset(subSpecific, GROUP==g & DMR==d)
#     totalPi<-rep(sum(subS$mean))
#     tmp<-cbind(subS, totalPi)
#     toPie<-rbind(toPie, tmp)
#   }
# }
# toPie$proportionPi<-toPie$mean*100/toPie$totalPi
# blank_theme <- theme_minimal()+
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))
# 
# ggplot(toPie, aes(x = "", y = proportionPi, fill = ANNOTATION)) +
#   geom_col() +
#   coord_polar(theta = "y")+
#   facet_grid(DMR~GROUP)+
#   scale_fill_manual(values = c("#ccbce2","#9d7ec7", "#d5d3d4","#b68090", "#f99299"))+
#   blank_theme+
#   theme(axis.text.x=element_blank())
# ggsave(paste0(outDir,"2_specific-epiPi_pie.png"),  width = 20, height = 20, units = "cm")

# or bar plot for the proportion
ggplot(toPie, aes(x=GROUP, y=proportionPi, fill=ANNOTATION))+
  geom_bar(stat='identity', position=position_dodge())+
  coord_flip()+ 
  facet_wrap(~DMR)+
  scale_fill_manual(values = c("#ccbce2","#9d7ec7", "#d5d3d4","#b68090", "#f99299"))+
  theme_minimal()
ggsave(paste0(outDir,"specific-epiPi_barProportion.png"),  width = 20, height = 20, units = "cm")
# or stacked bar plot for the proportion
ggplot(toPie, aes(x=GROUP, y=proportionPi, fill=ANNOTATION))+
  geom_bar(stat='identity')+
  coord_flip()+ 
  facet_wrap(~DMR)+
  scale_fill_manual(values = c("#ccbce2","#9d7ec7", "#d5d3d4","#b68090", "#f99299"))+
  theme_minimal()
ggsave(paste0(outDir,"specific-epiPi_barStackedProportion.png"),  width = 20, height = 20, units = "cm")

# or stacked bar plot for the mean PI
ggplot(toPie, aes(x=GROUP, y=mean, fill=ANNOTATION))+
  geom_bar(stat='identity')+
  coord_flip()+ 
  facet_wrap(~DMR)+
  scale_fill_manual(values = c("#ccbce2","#9d7ec7", "#d5d3d4","#b68090", "#f99299"))+
  theme_minimal()
ggsave(paste0(outDir,"specific-epiPi_barMeanPi.png"),  width = 20, height = 20, units = "cm")

######### along the chromosomes
global<-subset(epiPi, ANNOTATION=="general")
head(global)
global$BIN_START<-global$BIN_START/1000000 # in Mb
global$PIsmt<-global$PIsmt*1000 #10e3
dmr<-unique(global$DMR)
chr<-unique(global$CHROM)

subGlobal<-subset(global, GROUP!="3-SLC-wIntro" & 
                      GROUP!="3-SLC-woIntro" & 
                      GROUP!="4-SLL-wIntro" & 
                      GROUP!="4-SLL-woIntro")
for (d in dmr){
  tmp<-subset(subGlobal, DMR==d)
  for(c in chr){
    tmp2<-subset(tmp, CHROM==c)
    ggplot(tmp2, aes(x=BIN_START, y=PIsmt, group=GROUP)) +
      geom_line(aes(color=GROUP), linewidth=0.5)+
      facet_wrap(~CHROM)+
      #scale_x_continuous(limits = c(1, 20))+
      scale_color_manual(values=c("#209f68","#f6a29a","#d75e5e","#c33131"))
    ggsave(paste0(outDir,"3.1_smooth_",d,"_",c,"_epiAlongChr.png"),  width = 20, height = 20, units = "cm") 
  }
}

################ difference
general<-subset(epiPi, ANNOTATION=="general")
general$PI<-general$PI*1000 #10e3
general$BIN_START<-general$BIN_START/1000000 # in Mb

head(general)
gen<-general[,c(1,2,8,6,5)] # orden importa
head(gen)
tmp<-reshape2::dcast(gen, CHROM+BIN_START+DMR ~ GROUP, value.var = "PI")
head(tmp)
#dcast formula dcast(dataframe, IDvar1 +IDvar2 ~ ColVar, value.var= "colname")
# IDvar are those columns that will remain the same
# ColVar is the columns that will be the name of the new columns
# value.var will be the values of the nuevas columnas

ratio<-tmp[,c(1,2,3)]
ratio$WILD.PIM<-1-(round(tmp$`2-PIM`,3)/round(tmp$`1-WILD`,3)) # green to red transition
ratio$PIM.SLC<-1-(round(tmp$`3-SLC`,3)/round(tmp$`2-PIM`,3))
ratio$SLC.SLL<-1-(round(tmp$`4-SLL`,3)/round(tmp$`3-SLC`,3))
ratio<-reshape2::melt(ratio, id=c("CHROM", "BIN_START", "DMR"), variable.name="RATIO")
head(ratio)

chr<-unique(ratio$CHROM)
for (c in chr){
  toPlot<-subset(ratio, CHROM==c)
  xMax<-max(toPlot$BIN_START)
  toPlot['value'][toPlot['value']>1]<-1
  toPlot['value'][toPlot['value']<=-1]<-(-1)

  ggplot(toPlot, aes(x=BIN_START, y=value, fill=value<(-0.8))) +
    geom_bar(stat='identity')+
    scale_fill_manual(values=c("TRUE"="#ecae13", "FALSE"="grey"))+
    geom_line(aes(x = BIN_START, y = 0.8), linewidth = 0.5, color="red")+
    scale_x_continuous(breaks = seq(0,xMax, by=10))+
    scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1, by=0.2))+
    xlab(paste0("Position on chromosome ", gsub('SL2.50ch', "", c), " (Mb)"))+ ylab("Ratio")+
    facet_grid(DMR~RATIO)+
    theme_light()
  ggsave(paste0(outDir,"4a_allRatios_",c,".png"),  width = 20, height = 20, units = "cm")
}
