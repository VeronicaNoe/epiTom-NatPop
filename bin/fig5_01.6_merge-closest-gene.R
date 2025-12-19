suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
})
setwd("10_data-analysis/Fig5/ab_data-analysis/results")
outDir<-"10_data-analysis/Fig5/ab_data-analysis/results/"
input<-list.files(path=outDir, pattern = "closest", full.names = FALSE)
out<-c()
for(i in input){
  file_parts <- unlist(strsplit(i, "\\_"))  # Split on period
  ann <-gsub('-closest-','_',file_parts[2])
  dt<-data.table::fread(paste0(outDir,i),sep = '\t', data.table = T, header = F, 
                        fill = TRUE, na.string=c("NA"), nThread = 20)
  if(ann=="QTL"){
    annotation<-"QTL"
    dt$V12<-gsub('Name=', '', dt$V12)
    toKeep<-c(4,9,12,11)
  }else{
    annotation<-"QTLoverTE"
    dt$geneName<-annotation
    toKeep<-c(4,15,18,17)
  }
  dt<-dt[,..toKeep]
  colnames(dt)<-c('QTL', 'geneName','distance','geneFunction')
  dt$annotation<-annotation
  out<-rbind.data.frame(out,dt)
}
out
data.table::fwrite(out, file= paste0(outDir,'05.0_QTL-annotated.tsv'), quote=F,
                     row.names=F,col.names = T,sep="\t")
