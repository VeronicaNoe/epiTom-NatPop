library(dplyr)
library(ggplot2)
library(data.table)
#library(eulerr)
#library("RColorBrewer")
library(gplots)
library(ggplot2)
#library(viridis)
WD<-"/home/IPS2/vibanez/Desktop/Q-lab/2023/Fruits"
setwd(paste0(WD,'/data'))
outDir<-paste0(WD,"/results/")

df<- data.table::fread('leaves-fruit.meanMeth', sep = '\t', data.table = FALSE, fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
colnames(df)<-c('Sample','Tissue','ctxt','chr','meth','sd')
head(df)
toPlot<-subset(df, Tissue=="leaf" | Tissue=="fruit")
# box plot
  ggplot(toPlot, aes(x=Tissue, y=meth, fill=Tissue))+
    geom_boxplot()+
    scale_fill_manual(values=c("#c33131", "#209f68"))+
    facet_wrap(~ctxt) 
ggsave(paste0(outDir,"global-methylation_boxplot.png"),  width = 20, height = 20, units = "cm")
head(toPlot)
summ<-plyr::ddply(toPlot, c("Tissue","ctxt"), summarise,
                  mean= mean(meth))


100-((summ[3,'mean']*100)/summ[6,'mean'])
summ[3,'mean']
  summ[6,'mean']*2.52
####
sampleTab<-data.table::fread("/home/IPS2/vibanez/Desktop/Q-lab/2023/01_landscape/aa_data/00_sampleTab.tsv", header = TRUE,sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
head(sampleTab)
commonSample<-setdiff(sampleTab$Sample1, df$Sample)
############ DMR 
input<-list.files(pattern = "bed", full.names = FALSE)
df<- data.table::fread(input[14], sep = '\t', data.table = FALSE, fill = TRUE, check.names=FALSE,na.string=c("NA"), nThread = 10)
head(df)
