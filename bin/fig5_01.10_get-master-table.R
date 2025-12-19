suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

outDir<-"10_data-analysis/Fig5/ab_data-analysis/results/"
outPlot<-"10_data-analysis/Fig5/ab_data-analysis/plots/"

## load DMR annotation
mQTLanno<- data.table::fread(paste0(outDir,'05.0_QTL-annotated.tsv'),
                             sep = '\t', data.table = T, header = T,
                             fill = TRUE, na.string=c("NA"), nThread = 20)

mQTLanno$geneName <- sub("\\.\\d+$", "", mQTLanno$geneName)

## load meth-expression correlation
meth.exp<- data.table::fread(paste0(outDir,'06.2_all-correlation_expression-meth.tsv'),
                             sep = '\t', data.table = T, header = T,
                             fill = TRUE, na.string=c("NA"), nThread = 20)
head(meth.exp)
## load meth-metabolite correlation
meth.meta<- data.table::fread(paste0(outDir,'03.0_methylation-metabolite_correlation.tsv'),
                              sep = '\t', data.table = T, header = T,
                              fill = TRUE, na.string=c("NA"), nThread = 20)
head(meth.meta)
## load metabolite names
#
metaNames<- data.table::fread('01_raw-data/01.3_phenotypic-traits/metabolites_name.csv', sep = '\t', data.table = T, 
                              header = T,fill = TRUE, na.string=c("NA"), nThread = 20)
## load DMR target KO
koTargets<- data.table::fread(paste0(outDir,'07.1_DMRs-over-KO.tsv'), sep = '\t', 
                              data.table = T, header = T,fill = TRUE, na.string=c("NA"), nThread = 20)
metaNames[1:6,]
meth.exp[1:6,]
## merge metabolite names and mQTL
toKeep<-intersect(metaNames$Name, meth.exp$metabolite)
metaName_methExp <- metaNames %>%
  left_join(meth.exp, by = c("Name" = "metabolite")) %>%
  filter(Name %in% toKeep)

## merge with meth-exp correlation

metaName_methExp[1:6,1:6]
mQTLanno[1:6,1:4]

metaName_mQTLanno_corExpMeth <- metaName_methExp %>%
  left_join(mQTLanno, by = c("QTL" = "QTL","distance" = "distance","genes"="geneName", 
    "geneFunction" = "geneFunction"))

head(metaName_mQTLanno_corExpMeth)
## merge with meth-meta correlation
head(metaName_mQTLanno_corExpMeth)
meth.meta[1:6,]

metaName_mQTLanno_corExpMeth_corMethMeta <- metaName_mQTLanno_corExpMeth %>%
  left_join(meth.meta, by = c("Name" = "metabolite", "QTL" = "DMR"))

##  add ko targets
metaName_mQTLanno_corExpMeth_corMethMeta
koTargets


metaName_mQTLanno_corExpMeth_corMethMeta <- metaName_mQTLanno_corExpMeth %>%
  left_join(koTargets, by = c("QTL" = "DMR")) %>%
  select(-chr, -start,-end,-DMRtarget)

data.table::fwrite(metaName_mQTLanno_corExpMeth_corMethMeta, 
                   file= paste0(outDir,'08.1_master_table.tsv'), 
                   quote=F,row.names=F,col.names = T,sep="\t")
#### candidate genes
# close to the DMR
# expression affected
candidate_genes<-metaName_mQTLanno_corExpMeth_corMethMeta %>%
  filter(pValExpMeth_Pearson<=0.05 & abs(distance) <=10000 & kinship =="DMR-SNP")

data.table::fwrite(candidate_genes, file= paste0(outDir,'08.2_candidate_genes.tsv'), 
                   quote=F,row.names=F,col.names = T,sep="\t")

###############

