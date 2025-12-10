#!/users/bioinfo/vibanez/miniconda3/envs/methylKit/bin/Rscript
path2results<-"/mnt/data4/workspace3/vibanez/06_methylKitAnalysis/02_results/"
# get the folder passed from the shell
suppressPackageStartupMessages({
  library("methylKit")})
args <- commandArgs(trailingOnly = TRUE)
# check empty folders
allKOs<-list.files(path=args[1], pattern = "report.txt", full.names = FALSE)
if (length(allKOs) == 0) {
  stop("Empty folder. Skip.\n", call. = FALSE)
}

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (folder to process).\n", call. = FALSE)
} else {
  cat(paste0("====== Folder to process: ", args[1],"======"), sep="\n")
}


# list all files
chr2process=unlist(strsplit(allKOs[1],"_"))
chr2process=grep('chr', chr2process, value=TRUE)
chr2process=sub(".report.txt", "",chr2process)

allWT<-list.files(path= '/mnt/data4/workspace3/vibanez/06_methylKitAnalysis/01_rawData/WT/',pattern = "Ohio", full.names = FALSE)
WT<-grep(chr2process, allWT, value = TRUE)

# print info in screen
cat("----- Input files:", sep="\n")
cat(allKOs[1], sep="\n")
cat("...", sep="\n")
cat(allKOs[length(allKOs)], sep="\n")
cat(WT, sep="\n")
accession<-c(WT,allKOs)
#add path to Control samples
accession<-sub('Ohio', '/mnt/data4/workspace3/vibanez/06_methylKitAnalysis/01_rawData/WT/Ohio', accession)
#cat(accession, sep="\n")

inControl<-grep("_CHG_", WT, value = TRUE)
inControl<-sub(paste0("_CHG_",chr2process,".report.txt"), "",inControl)
inKO<-grep("_CHG_", allKOs, value = TRUE)
inKO<-sub(paste0("_CHG_",chr2process,".report.txt"), "",inKO)
#add path to KO samples
accession<-sub(inKO, paste0('/mnt/data4/workspace3/vibanez/06_methylKitAnalysis/01_rawData/KO/',inKO,"/",chr2process,"/",inKO), accession)

allSamples<-c(inControl,inKO)  

sampleNames<-as.list(c(allSamples)) #Definir el nombre de las muestras para methylkit
  
ctrlContrast<-c()
for(i in 1:length(inControl)){
  ctrlContrast<-append(ctrlContrast,0)
}
sampleContrast<-c()
for(i in 1:length(inKO)){
  sampleContrast<-append(sampleContrast, 1)
}
contr<-c(ctrlContrast, sampleContrast) #Definir los contrastas del diseño exp
refGen<-"SL2.50" #Nombre del genoma de Referencia
  
CGdata<-list()
CHGdata<-list()
CHHdata<-list()
inFileCG<-grep("CG", c(accession), value=TRUE)
inFileCHG<-grep("CHG", c(accession), value=TRUE)
inFileCHH<-grep("CHH", c(accession), value=TRUE)

for (i in 1:length(inFileCG)){
    CGdata[[i]] <- inFileCG[i]
    CHGdata[[i]] <- inFileCHG[i]
    CHHdata[[i]] <- inFileCHH[i]
}

#start with methylkit
CG.cov<-methRead(CGdata, context = "CpG", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                   mincov = 0, header = FALSE)
  
CHG.cov<-methRead(CHGdata, context = "CHG", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                  mincov = 0, header = FALSE)
  
CHH.cov<-methRead(CHHdata, context = "CHH", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                  mincov = 0, header = FALSE)

CG.filtered=filterByCoverage(CG.cov,lo.count=2,lo.perc=NULL,
                               hi.count=NULL,hi.perc=99.9)
  
CHG.filtered=filterByCoverage(CHG.cov,lo.count=2,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)
  
CHH.filtered=filterByCoverage(CHH.cov,lo.count=2,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)
  
CG.norm<-unite(CG.filtered, destrand = FALSE, min.per.group = 1L)
CHG.norm<-unite(CHG.filtered, destrand = FALSE, min.per.group = 1L)
CHH.norm<-unite(CHH.filtered, destrand = FALSE, min.per.group = 1L)
  
## Methylation calling en bins (o tiles) en lugar de por base para cada (se puede variar el win size)
CG.wind <- tileMethylCounts(CG.norm, win.size = 100, step.size = 100)
CHG.wind <- tileMethylCounts(CHG.norm, win.size = 100, step.size = 100)
CHH.wind <- tileMethylCounts(CHH.norm, win.size = 100, step.size = 100)
  
CG.bins <- calculateDiffMeth(CG.wind, mc.cores= 20)
CHG.bins <- calculateDiffMeth(CHG.wind, mc.cores= 20)
CHH.bins <- calculateDiffMeth(CHH.wind, mc.cores = 28)
  
# DMR with 40% of difference and qvalue 0.05
CG.DMR_40_0.05.hyper=getMethylDiff(CG.bins, difference=40, qvalue=0.05,type="hyper")
CG.DMR_40_0.05.hypo=getMethylDiff(CG.bins, difference=40, qvalue=0.05,type="hypo")
  
CHG.DMR_40_0.05.hyper=getMethylDiff(CHG.bins, difference=20, qvalue=0.05,type="hyper")
CHG.DMR_40_0.05.hypo=getMethylDiff(CHG.bins, difference=20, qvalue=0.05,type="hypo")
  
CHH.DMR_40_0.05.hyper=getMethylDiff(CHH.bins, difference=10, qvalue=0.05,type="hyper")
CHH.DMR_40_0.05.hypo=getMethylDiff(CHH.bins, difference=10, qvalue=0.05,type="hypo")
  
# save files  
write.table(data.frame(CG.DMR_40_0.05.hyper), file=paste0(path2results,"DMR_hyper_CG_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CG.DMR_40_0.05.hypo), file=paste0(path2results,"DMR_hypo_CG_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHG.DMR_40_0.05.hyper), file=paste0(path2results,"DMR_hyper_CHG_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHG.DMR_40_0.05.hypo), file=paste0(path2results,"DMR_hypo_CHG_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHH.DMR_40_0.05.hyper), file=paste0(path2results,"DMR_hyper_CHH_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHH.DMR_40_0.05.hypo), file=paste0(path2results,"DMR_hypo_CHH_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
# save coverage data after filtering
write.table(data.frame(CG.wind), file=paste0(path2results,"00_filterd_CG_wind_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHG.wind), file=paste0(path2results,"00_filterd_CHG_wind_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHH.wind), file=paste0(path2results,"00_filterd_CHH_wind_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
#save bin files with q values
write.table(data.frame(CG.bins), file=paste0(path2results,"00_filterd_CG_bins_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHG.bins), file=paste0(path2results,"00_filterd_CHG_bins_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(data.frame(CHH.bins), file=paste0(path2results,"00_filterd_CHH_bins_",inKO,"_",chr2process,".txt"), quote=F,row.names=F,sep="\t")
write.table(sampleNames, file=paste0(inKO,"_",chr2process,"_sampleNames.tsv"), quote=F,row.names=F,sep="\t")
