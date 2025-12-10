#!/users/bioinfo/vibanez/miniconda3/envs/methylKit/bin/Rscript
path2results<-"/mnt/data4/workspace3/vibanez/01_WT-KO-analysis/02_description_WT-KO/results/"


# get the folder passed from the shell
suppressPackageStartupMessages({
  library("methylKit")})
args <- commandArgs(trailingOnly = TRUE)
# check empty folders
allSamples<-list.files(path=args[1], pattern = "report.txt", full.names = FALSE)
#allSamples<-list.files(pattern = "report.txt", full.names = FALSE)
ohio<-grep("Ohio", allSamples, value=TRUE, invert = FALSE)
met<-grep("met1", allSamples, value=TRUE, invert = FALSE)
allSamples<-c(ohio, met)
cat(allSamples)
# chromosome to process
chr2process=unlist(strsplit(allSamples[1],"_"))
chr2process=grep('chr', chr2process, value=TRUE)
chr2process=sub(".report.txt", "",chr2process)
path2sample<-paste0("/mnt/data4/workspace3/vibanez/01_WT-KO-analysis/02_description_WT-KO/data/",chr2process,"/",allSamples)


# define names for each element in the list
accessions<-grep("_CHG", allSamples, value=TRUE)
accessions<-sub(paste0("_CHG_",chr2process,".report.txt"), "",accessions)
sampleNames<-as.list(c(accessions)) #Definir el nombre de las muestras

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
CG.cov<-methRead(CGdata, context = "CpG", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                   mincov = 0, header = FALSE)
  
CHG.cov<-methRead(CHGdata, context = "CHG", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                  mincov = 0, header = FALSE)
  
CHH.cov<-methRead(CHHdata, context = "CHH", sample.id =sampleNames,
                  treatment =contr, assembly= refGen , pipeline = "bismarkCytosineReport",
                  mincov = 0, header = FALSE)

# make the table
CG.df<-data.frame(unite(CG.cov, destrand = FALSE, min.per.group = 0L))
ctxt<-rep("CG", times=nrow(CG.df))
CG<-cbind(CG.df, ctxt)

CHG.df<-data.frame(unite(CHG.cov, destrand = FALSE, min.per.group = 0L))
ctxt<-rep("CHG", times=nrow(CHG.df))
CHG<-cbind(CHG.df, ctxt)

CHH.df<-data.frame(unite(CHH.cov, destrand = FALSE, min.per.group = 0L))
ctxt<-rep("CHH", times=nrow(CHH.df))
CHH<-cbind(CHH.df, ctxt)

# unite all
out<-rbind(CG, CHG, CHH)
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
colnames(out)<-toName

# save files  
data.table::fwrite(out, file=paste0(path2results, chr2process,"_Ohio-WT-KO.methylation.bed"), quote=F,row.names=F,col.names = T,sep="\t")
