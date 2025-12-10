#!/usr/bin/R
suppressPackageStartupMessages({
	require("data.table")
	library(ggplot2)
	library(dplyr)})
#install.packages("dplyr")
#library(dplyr)
#color for tomato groups
args <- commandArgs(trailingOnly = TRUE)
wd<-"/mnt/disk2/vibanez/02_methylkit/ae_cumulative-DMR"
input<-paste0(wd, "/bb_output/")
output<-paste0(wd,"/bc_results/")
infile<-paste0(input,args[1])
#print(infile)

df<- data.table::fread(infile, sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
infoCol<-c("chr", "start", "end", "V133")
toKeep<-setdiff(colnames(df), infoCol)
#print(toKeep)

df<-df[,toKeep]
#colnames(df)
#write.table(df, paste0(output,args[1],"_check"),quote=F,row.names=F,sep="\t")

## start the loop 
nAcc=ncol(df)
nPermutation=100
outDF<-c()

for(a in 1:nAcc){
	print(a)
	for(p in 1:nPermutation){
		if(a==1){
			toKeep<-sample(1:ncol(df),a,  replace = FALSE)
			tmp<-as.data.frame(df[,toKeep])
			tmp$toRemove<-is.na(tmp)
			tmp<-subset(tmp, toRemove!='TRUE')
			nDMR<-nrow(tmp)
			out<-cbind.data.frame(a,p,nDMR)
			outDF<-rbind.data.frame(outDF,out)
		} else {
			toKeep<-sample(1:ncol(df),a,  replace = FALSE)
			tmp<-df[,toKeep]
			tmp$toRemove<-rowSums(is.na(tmp))==a
			tmp<-subset(tmp, toRemove!='TRUE')
			nDMR<-nrow(tmp)
			out<-cbind.data.frame(a,p,nDMR)
			outDF<-rbind.data.frame(outDF,out)
		}
	}
}
write.table(outDF, paste0(output,args[1],"_sampling.tsv"),quote=F,row.names=F,sep="\t")
summ<-plyr::ddply(outDF, c("a"), summarise,
	mean = mean(nDMR),
	sd = sd(nDMR),
	sem = sd(nDMR)/sqrt(length(nDMR)))
write.table(summ, paste0(output, args[1],"_summary-sampling.tsv"),quote=F,row.names=F,sep="\t")
#pdf(paste0(output, args[1],"_sampling.pdf"))
#ggplot(summ, aes(x = a, y = mean)) +
#  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), 
#              alpha = 0.5, color = "grey") +
#  geom_line()
#dev.off()
