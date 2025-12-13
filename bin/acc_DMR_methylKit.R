#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library("methylKit")})
args <- commandArgs(trailingOnly = TRUE)
/03_biseq-processing/03.4_filtering/ac_filter
path2acc<-"/03_biseq-processing/03.4_filtering/ac_filter"
path2results<-"04_methylome-comparison/aa_natural-accessions"
chr2process=args[2]

#
acc<-list.files(path=path2acc, pattern =paste0(chr2process,".filtered.bed") , full.names = FALSE)
acc<-grep(paste0(args[1],"_Biseq_C"), acc, value=TRUE)
#acc<-grep("LA0534", acc, value=TRUE)
acc<-grep("filterOut", acc,invert=TRUE, value=TRUE)

ref2compare<-list.files(path= path2acc, pattern = "TS-253_leaf", full.names = FALSE)
ref2compare<-grep(paste0(chr2process,".filtered.bed"), ref2compare, value =TRUE )
ref2compare<-grep("filterOut", ref2compare, invert=TRUE, value=TRUE)
allacc<-c(ref2compare, acc)
allacc<-paste0(path2acc, allacc)
accName<-args[1]
allacc

allSamples<-grep("CG",allacc, value = TRUE)
allSamples<-gsub(path2acc, "", allSamples)
sampleNames<-as.list(c(allSamples)) #Definir el nombre de las muestras para methylkit
ctrlContrast<-0
sampleContrast<-1
contr<-c(ctrlContrast, sampleContrast) #Definir los contrastas del diseño exp
refGen<-"SL2.50" #Nombre del genoma de Referencia

CGdata<-list()
CHGdata<-list()
CHHdata<-list()
inFileCG<-grep("CG", c(allacc), value=TRUE)
inFileCHG<-grep("CHG", c(allacc), value=TRUE)
inFileCHH<-grep("CHH", c(allacc), value=TRUE)

for (i in 1:length(inFileCG)){
  CGdata[[i]] <- inFileCG[i]
  CHGdata[[i]] <- inFileCHG[i]
  CHHdata[[i]] <- inFileCHH[i]
}

#start with methylkit
cat(paste0("----- Loading data to mKit", "\n"))
cat(paste0("---------- CG", "\n"))
CG.cov<-methRead(CGdata, context = "CpG", sample.id =sampleNames,
                 treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                 mincov = 0, header = FALSE)
cat(paste0("---------- CHG", "\n"))
CHG.cov<-methRead(CHGdata, context = "CHG", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                  mincov = 0, header = FALSE)
cat(paste0("---------- CHH", "\n"))
CHH.cov<-methRead(CHHdata, context = "CHH", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                  mincov = 0, header = FALSE)

cat(paste0("----- Filtering high coverage", "\n"))
cat(paste0("---------- CG", "\n"))
CG.filtered=filterByCoverage(CG.cov,lo.count=5,lo.perc=NULL,
                             hi.count=NULL,hi.perc=99.9, chunk.size = 1e+10)
cat(paste0("---------- CHG", "\n"))
CHG.filtered=filterByCoverage(CHG.cov,lo.count=5,lo.perc=NULL,
                              hi.count=NULL,hi.perc=99.9, chunk.size = 1e+16)
cat(paste0("---------- CHH", "\n"))
CHH.filtered=filterByCoverage(CHH.cov,lo.count=5,lo.perc=NULL,
                              hi.count=NULL,hi.perc=99.9, chunk.size = 1e+20)

cat(paste0("----- Merging data", "\n"))
cat(paste0("---------- CG", "\n"))
CG.norm<-unite(CG.filtered, destrand = FALSE, min.per.group = 1L)
cat(paste0("---------- CHG", "\n"))
CHG.norm<-unite(CHG.filtered, destrand = FALSE, min.per.group = 1L)
cat(paste0("---------- CHH", "\n"))
CHH.norm<-unite(CHH.filtered, destrand = FALSE, min.per.group = 1L)

cat(paste0("----- Creating windows", "\n"))
cat(paste0("---------- CG", "\n"))
CG.wind <- tileMethylCounts(CG.norm, win.size = 100, step.size = 100, mc.cores = 50)
cat(paste0("---------- CHG", "\n"))
CHG.wind <-tileMethylCounts(CHG.norm, win.size = 100, step.size = 100, mc.cores= 50)
cat(paste0("---------- CHH", "\n"))
CHH.wind <- tileMethylCounts(CHH.norm, win.size = 100, step.size = 100, mc.cores = 50)

cat(paste0("----- Calculating DMR", "\n"))
cat(paste0("---------- CG", "\n"))
CG.DMR <- data.frame(calculateDiffMeth(CG.wind, mc.cores= 20))
cat(paste0("---------- CHG", "\n"))
CHG.DMR <- data.frame(calculateDiffMeth(CHG.wind, mc.cores= 20))
cat(paste0("---------- CHH", "\n"))
CHH.DMR <- data.frame(calculateDiffMeth(CHH.wind, mc.cores = 28))

control<-unlist(strsplit(ref2compare[1],"_"))[1]
toNameCol<-c(control,accName)
cnames<-c()
for (i in toNameCol){
  total<-paste0(i,'_total')
  meth<-paste0(i,'_methylated')
  unmeth<-paste0(i,'_unmethylated')
  tmp<-c(total, meth, unmeth)
  cnames<-c(cnames, tmp)
}

cat(paste0("----- Filtering by qvalue and methylation difference", "\n"))
cat(paste0("---------- CG", "\n"))
CG.dfDMR<-data.frame(CG.wind)
context<-rep("CG", times=nrow(CG.dfDMR))
toName<-c('chr', 'start', 'end', 'strand', cnames)
colnames(CG.dfDMR)<-toName
CG_DMR<-cbind.data.frame(CG.dfDMR, context, CG.DMR[,5:7])
CG_DMR.toSave<-subset(CG_DMR, qvalue<=0.05 & (meth.diff <= -50 | meth.diff >= 50))
data.table::fwrite(CG_DMR.toSave, file=paste0(path2results, accName,"_",chr2process,"_CG_DMR.methylation.bed"),
                   quote=F,row.names=F,col.names = F,sep="\t", nThread=20, buffMB=1024)

cat(paste0("---------- CHG", "\n"))
CHG.dfDMR<-data.frame(CHG.wind)
context<-rep("CHG", times=nrow(CHG.dfDMR))
toName<-c('chr', 'start', 'end', 'strand', cnames)
colnames(CHG.dfDMR)<-toName
CHG_DMR<-cbind.data.frame(CHG.dfDMR, context, CHG.DMR[,5:7])
CHG_DMR.toSave<-subset(CHG_DMR, qvalue<=0.05 & (meth.diff <= -50 | meth.diff >= 50))
data.table::fwrite(CHG_DMR.toSave, file=paste0(path2results, accName,"_",chr2process,"_CHG_DMR.methylation.bed"),
                   quote=F,row.names=F,col.names = F,sep="\t", nThread=20, buffMB=1024)

cat(paste0("---------- CHH", "\n"))
CHH.dfDMR<-data.frame(CHH.wind)
context<-rep("CHH", times=nrow(CHH.dfDMR))
toName<-c('chr', 'start', 'end', 'strand', cnames)
colnames(CHH.dfDMR)<-toName
CHH_DMR<-cbind.data.frame(CHH.dfDMR, context, CHH.DMR[,5:7])
CHH_DMR.toSave<-subset(CHH_DMR, qvalue<=0.05 & (meth.diff <= -10 | meth.diff >= 10))
data.table::fwrite(CHH_DMR.toSave, file=paste0(path2results, accName,"_",chr2process,"_CHH_DMR.methylation.bed"),
                   quote=F,row.names=F,col.names = F,sep="\t", nThread=20, buffMB=1024)

