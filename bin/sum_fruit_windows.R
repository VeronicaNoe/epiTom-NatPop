#!/users/bioinfo/vibanez/miniconda3/envs/methylKit/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
path2acc<-"08_fruit-processing/08.1_DMR-classification/08.1_get-collapsed-loci/tmp/"
path2results<-"08_fruit-processing/08.1_DMR-classification/08.1_get-collapsed-loci"
#args<-"ch01_S-lycLYC1410.C-DMR.bed"
data<- data.table::fread(paste0(path2acc, args[1]), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20) 
toKeep<-unique(data$V2)
outFile<-c()
for (i in toKeep){
	subdf<-subset(data, V2==i)
	if(dim(subdf)[1]!=1){
		same<-c("V1", "V2", "V3", "V4")
		V8=sum(subdf[1:nrow(subdf),"V8"])
		V9=sum(subdf[1:nrow(subdf),"V9"])
		V10=sum(subdf[1:nrow(subdf),"V10"])
		#V12=mean(subdf[1:nrow(subdf),"V12"])
		#V13=mean(subdf[1:nrow(subdf),"V13"])
		#V14=mean(subdf[1:nrow(subdf),"V14"])
		 V11=paste0(subdf[1,"V11"],"-", subdf[2,"V11"],"-", subdf[3,"V11"])
		#tmp<-cbind.data.frame(subdf[1,c(same)],V8,V9,V10,V11,V12,V13,V14)
		tmp<-cbind.data.frame(subdf[1,c(same)],V8,V9,V10,V11)
		outFile<-rbind.data.frame(outFile, tmp)
	} else {
		outFile<-rbind.data.frame(outFile, subdf[,c(1,2,3,4,8,9,10,11)])
		#outFile<-rbind.data.frame(outFile, subdf[,c(1,2,3,4,8,9,10,11,12,13,14)])

	}
}
data.table::fwrite(outFile, file=paste0(path2results, args[1]), quote=F,row.names=F,col.names = F,sep="\t", nThread=20, buffMB=1024)
