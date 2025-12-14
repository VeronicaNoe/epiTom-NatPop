inDir<-"06_get-meth-vcf/ab_epialleles/"
outDir<-("06_get-meth-vcf/ac_vcf-metabolome/")

dmrs<-c("C-DMR","CG-DMR")
out<-c()
for( d in dmrs ){
  df<- data.table::fread(paste0(inDir,d,"_general_leaf-metabolome_LD.vcf"), 
                         sep = '\t', data.table = F, skip="#CHROM",
                         fill = TRUE, na.string="NA", nThread = 20)
  colnames(df)<-gsub("0_","", colnames(df))
  out<-rbind.data.frame(out, df)
}
out <- out[order(out$`#CHROM`, as.numeric(out$POS)), ]
data.table::fwrite(out, file= paste0(outDir,'DMR_general_leaf-metabolome_LD.vcf'),
                   quote=F,row.names=F,col.names = T,sep="\t")

#######
#LD file
snps <- data.table::fread("",
                          skip = "#CHROM", sep = "\t", data.table = FALSE, fill = TRUE, header = TRUE)

out$`#CHROM` <- as.integer(gsub('SL2.50ch','',out$`#CHROM`))
toKeep<-intersect(colnames(snps),colnames(out))
out<-out[,toKeep]
snps<-snps[,toKeep]
# merge
dmr_snp<-rbind.data.frame(snps, out)
data.table::fwrite(dmr_snp, file= paste0(outDir,'DMR-SNP_general_leaf-metabolome_LD.vcf'), 
                   quote=F, row.names=F, col.names = T, sep="\t")
