## run this in mclintock
setwd("/mnt/disk2/vibanez/metabolome-GWAS/leaf/checkGIF")
indir<-"/mnt/disk2/vibanez/metabolome-GWAS/leaf/checkGIF/data/"
outPlot<-"/mnt/disk2/vibanez/metabolome-GWAS/leaf/checkGIF/plots/"
library(qqman)
args <- commandArgs(trailingOnly = TRUE)
#input<-list.files(pattern = ".aBN.ps.gz", full.names = FALSE)
metabolite<-unlist(strsplit(args[1], '.', fixed=TRUE))[1]
molmarker<-unlist(strsplit(args[1], '.', fixed=TRUE))[3]
toPlot<-data.table::fread(paste0(indir, args[1]), sep = '\t',
                           data.table = FALSE, fill = T,header = F, 
                           na.string=c("NA"), nThread = 20)
  #format for the library
toPlot$V1<-gsub('SL2.50ch','', toPlot$V1)
toPlot$CHR <- sapply(strsplit(as.character(toPlot$V1), ":"),
                     function(x) as.integer(x[1]))
toPlot$BP <- sapply(strsplit(as.character(toPlot$V1), ":"),
                    function(x) as.integer(x[2]))
toPlot<-toPlot[,c(1,5,6,4)]
colnames(toPlot)<-c("SNP","CHR", "BP", "P")
nAss<-nrow(subset(toPlot, P<=5e-8))
# get lambda genomic inflation
chisq <-qchisq(1-toPlot$P, 1)
lambda<-median(chisq)/qchisq(0.5,1)
gif<-rbind.data.frame(metabolite, molmarker,nAss,lambda)

inflationFactor <- paste("Inflation Factor:", round(lambda, 5))
toPlot$CHR<-as.numeric(toPlot$CHR)
toPlot$BP<-as.numeric(toPlot$BP)
# make the plot
png(paste0(outPlot, metabolite,"_manhattanPlot_", molmarker, "_NonSig-by-lambda", ".png"))
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
manhattan(toPlot)
qq(toPlot$P)
middle_x <- par("pin")[1] - par("pin")[1] # width
middle_y <- (par("pin")[2] + par("pin")[2])*2  # height
text(x = middle_x, y = middle_y, label = inflationFactor, pos = 4,
     col = "red", cex = 1.0)
dev.off()

data.table::fwrite(gif, file= paste0(metabolite,'_',molmarker,'_lambda.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")
