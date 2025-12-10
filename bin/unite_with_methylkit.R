#!/users/bioinfo/vibanez/miniconda3/envs/methylKit/bin/Rscript
path2acc<-"/mnt/disk2/vibanez/02_methylkit/ab_mergeDMRs/aa_data/"
path2results<-"/mnt/disk2/vibanez/02_methylkit/ab_mergeDMRs/bb_output/"
# get the folder passed from the shell
suppressPackageStartupMessages({
  library("methylKit")})
args <- commandArgs(trailingOnly = TRUE)
# check empty folders
allSamples<-list.files(path=path2acc, pattern = args[1], full.names = FALSE)
#cat(paste0(allSamples, "\n"))

# chromosome to process
chr2process=args[1]
path2sample<-paste0(path2acc,allSamples)

# define names for each element in the list
accessions<-grep("_CHG", allSamples, value=TRUE)
sampleNames<-as.list(c(accessions)) #Definir el nombre de las muestras
#cat(paste0(accessions, "\n"))

# define contrasts
ctrlContrast<-rep(0, times=length(sampleNames)-1)
sampleContrast<-1
contr<-c(ctrlContrast, sampleContrast) #Definir los contrastas del diseño exp

#define reference genome
refGen<-"SL2.50" #Nombre del genoma de Referencia

# defin samples for each context
CGdata<-list()
CHGdata<-list()
CHHdata<-list()
inFileCG<-grep("CG", c(path2sample), value=TRUE)
inFileCHG<-grep("CHG", c(path2sample), value=TRUE)
inFileCHH<-grep("CHH", c(path2sample), value=TRUE)

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

#add col names
cnames<-c()
for (i in accessions){
  total<-paste0(i,'_total')
  meth<-paste0(i,'_methylated')
  unmeth<-paste0(i,'_unmethylated')
  tmp<-c(total, meth, unmeth)
  cnames<-c(cnames, tmp)
}
toName<-c('chr', 'start', 'end', 'strand', cnames, 'context')

# make the table
cat(paste0("----- Windows", "\n"))
cat(paste0("---------- CG", "\n"))
CG.wind <- tileMethylCounts(CG.cov, win.size = 100, step.size = 100, mc.cores = 50)
cat(paste0("---------- CHG", "\n"))
CHG.wind <-tileMethylCounts(CHG.cov, win.size = 100, step.size = 100, mc.cores= 50)
cat(paste0("---------- CHH", "\n"))
CHH.wind <- tileMethylCounts(CHH.cov, win.size = 100, step.size = 100, mc.cores = 50)

cat(paste0("----- Uniting", "\n"))
cat(paste0("---------- CG", "\n"))
CG.df<-data.frame(unite(CG.wind, destrand = FALSE, min.per.group = 0L))
ctxt<-rep("CG", times=nrow(CG.df))
CG<-cbind(CG.df, ctxt)
colnames(CG)<-toName
data.table::fwrite(CG, file=paste0(path2results, chr2process,""), quote=F,row.names=F,col.names = T,sep="\t")

cat(paste0("---------- CHG", "\n"))
CHG.df<-data.frame(unite(CHG.wind, destrand = FALSE, min.per.group = 0L))
ctxt<-rep("CHG", times=nrow(CHG.df))
CHG<-cbind(CHG.df, ctxt)
colnames(CHG)<-toName
data.table::fwrite(CHG, file=paste0(path2results, chr2process,""), quote=F,row.names=F,col.names = T,sep="\t")

cat(paste0("---------- CHH", "\n"))
CHH.df<-data.frame(unite(CHH.wind, destrand = FALSE, min.per.group = 0L))
ctxt<-rep("CHH", times=nrow(CHH.df))
CHH<-cbind(CHH.df, ctxt)
colnames(CHH)<-toName
data.table::fwrite(CHH, file=paste0(path2results, chr2process,""), quote=F,row.names=F,col.names = T,sep="\t")

