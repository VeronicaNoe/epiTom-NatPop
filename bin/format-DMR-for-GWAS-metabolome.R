# ACTIVATE WORKING
setwd("/mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles")
# get sample Names
sampleNames<-data.table::fread("/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype/raw/samples2get",
                               sep = '\t', data.table = FALSE, fill = TRUE,
                               header = T, na.string=c("NA"), nThread = 20)
rownames(sampleNames)<-sampleNames$forMeth
infoCol<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
# get input files
input<-list.files(pattern = "gz", full.names = FALSE)
#input<-grep('tbi',input, value = T, invert = T)

for( i in input ){
  dmr<-unlist(strsplit(i, '.', fixed=TRUE))[1]
  df<- data.table::fread(i, sep = '\t', data.table = F, skip="#CHROM",
                         fill = TRUE, na.string="NA", nThread = 20) 
  colnames(df)<-gsub('_leaf','',colnames(df))
  colnames(df)<-gsub('LA2838A','S-lycLA2838A', colnames(df))
  df$ID <- paste0(gsub('SL2.50','',df$ID),":",dmr)
  methName<-setdiff(colnames(df),infoCol)
  toKeep<-intersect(sampleNames$forMeth, methName)
  toRename<-sampleNames[toKeep,]
  df <- df[,c(infoCol,toKeep)]
  df$`#CHROM` <- as.integer(gsub('SL2.50ch','',df$`#CHROM`))
  colnames(df)<-c(infoCol,toRename$forSnps)
  df <- df[order(df$`#CHROM`, as.numeric(df$POS)), ]
  data.table::fwrite(df, file= paste0(dmr,'_general_leaf-metabolome.vcf'),
                     quote=F,row.names=F,col.names = T,sep="\t")
}