#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/bb_get-meth-levels-around-gene-windows")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/bc_format-data/"

# get colnames form RNAseq data
rnaData<-data.table::fread("/mnt/disk2/vibanez/07_rnaseq-processing/07.3_counts/leaf_transcriptome-counts_normalized_batch-adjusted.tsv", 
                           sep = '\t',
                           data.table = F, fill = TRUE, check.names=FALSE,
                           na.string=c("NA"), nThread = 10)
#
sampleNames<-data.table::fread("/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation/C-DMR_colNames.tsv",
                               sep = '\t',header=F,data.table = FALSE, fill = TRUE, check.names=FALSE,
                               na.string=c("NA"), nThread = 10)

sampleCol<-gsub('_leaf','', sampleNames$V1)
infoCol<-c('toRemove1','toRemove2','toRemove3','strand','geneName')
toNameCol<-c(infoCol,'toRemove4','chr','start','end',sampleCol)
toKeep<-intersect(colnames(rnaData), sampleCol)
rm(rnaData)
# load data
input<-list.files(pattern = "DMR", full.names = FALSE)
methLevels<-c()
for (i in input){
  annotation<-unlist(strsplit(i, '.', fixed=TRUE))[1]
  dmr<-unlist(strsplit(i, '.', fixed=TRUE))[2]
  df<-data.table::fread(i, sep = '\t',
                        data.table = F, fill = TRUE, check.names=FALSE,
                        na.string=c("NA"), nThread = 10)
  colnames(df)<-toNameCol
  df$dmr<-dmr
  df$genePart<-annotation
  df<-df[,c('chr','start','strand','geneName','genePart','dmr',toKeep)]
  methLevels<-rbind.data.frame(methLevels, df)
}
methLevels$toFilter<-paste0(gsub('SL2.50','',methLevels$chr),'_',methLevels$dmr,"_",methLevels$start)
methLevels$dmr<-NULL

methLevels<-reshape2::melt(methLevels, id=c('chr','start','toFilter','geneName','strand','genePart'),
                        variable.name="Sample", value.name = "methylationLevels")
#head(methLevels)

############################################################
# split methylation levels by chr and dmr
############################################################
byChr<-unique(methLevels$chr)
dmr<-c("C-DMR","CG-DMR")
for( c in byChr){
  tmp<-subset(methLevels, chr == c)
  for( d in dmr){
    out<-tmp[grepl(d, tmp$toFilter),]
    data.table::fwrite(out, 
                       file= paste0(outDir,gsub('SL2.50','',c),"_",d,'_methylation-levels.tsv'), quote=F,
                       row.names=F,col.names = T,sep="\t")
  }
}
