basedir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype"
setwd(paste0(basedir,"/raw"))
outdir<-paste0(basedir,"/formated-files/")
df<- data.table::fread("leaf.metabolites.tsv", sep = '\t', data.table = FALSE,
        fill = TRUE,na.string="NA", nThread = 20)

toKeep<-grep('Accessions', colnames(df), invert = T, value = T)
df<-cbind.data.frame(df$Accessions,log2(df[,toKeep]))
rownames(df)<-df[,1]
sampleNames<-data.table::fread("/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/SNP_general_leaf_LD.tfam",
                               sep = ' ', data.table = FALSE, fill = TRUE,
                               header = FALSE, na.string=c("NA"), nThread = 20)
df_sored<-df[sampleNames$V2,]
for(c in 2:ncol(df)){
        tmp<-df_sored[,c(1,1,c)]
        name2save<-colnames(df_sored)[c]
        data.table::fwrite(tmp, file=paste0(outdir, name2save,".leaf-metabolites"),
        na= "NA", quote=F,row.names=F,col.names = F,sep="\t")
}
