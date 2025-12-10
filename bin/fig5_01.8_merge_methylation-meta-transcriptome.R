suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(plotly)
  library(ggfortify)
  library(gridExtra)
})

outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/ab_data-analysis/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/ab_data-analysis/plots/"

  ############################################################
  ## load genes  & mQTL
  ############################################################
transcriptome<- data.table::fread("/mnt/disk2/vibanez/07_rnaseq-processing/07.3_counts/leaf-transcriptome_deseq-combat.tsv", 
                                                 sep = '\t', data.table = T, fill = TRUE, 
                                                 na.string=c("NA"), nThread = 20)


QTLannotated<- data.table::fread(paste0(outDir,"05.0_QTL-annotated.tsv"), 
                                    sep = '\t', data.table = T, fill = TRUE, 
                                    na.string=c("NA"), nThread = 20)
QTLannotated$geneName <- gsub("Name=", "", QTLannotated$geneName)
QTLannotated$geneName <- sub("\\.\\d+$", "", QTLannotated$geneName)
QTLannotated
  
transcriptome<-transcriptome[transcriptome$Geneid %in% QTLannotated$geneName,]
setDT(transcriptome)
geneEpx <- melt(transcriptome, id.vars = "Geneid", variable.name = "Accessions", value.name = "Expression")

merged_data <- QTLannotated %>%
    left_join(geneEpx, by = c("geneName" = "Geneid"))
  ############################################################
  ## load methylation levles
  ############################################################
meth_meta<-data.table::fread(paste0(outDir,'02.0_merged-mLevels-mStatus-metaLevels.tsv'), 
                               sep = '\t', data.table = T,
                               fill = TRUE, na.string=c("NA"), nThread = 20)
head(meth_meta)
  
merged_data <- merged_data %>%
  left_join(meth_meta, by = c("QTL" = "ID", "Accessions" = "Accessions")) %>%
    arrange()
  
  merged_data$expressionLog2<-log2(merged_data$Expression+1)
  toSort<-c('QTL','mm','kinship','resultType', 'annotation','metabolite','Accessions',
            'methStatus', 'methLevels', 'metaLevel','geneName',
            'Expression','expressionLog2','distance','geneFunction')
  
  merged_data<-merged_data[,..toSort]
  
  data.table::fwrite(merged_data, 
                     file= paste0(outDir,'05.2_expression-meth-meta.tsv'), quote=F,
                     row.names=F,col.names = T,sep="\t")
## plots
  genes<-unique(merged_data$geneName)
  for(g in 1:length(genes)){
    tmp<-subset(merged_data, geneName==genes[g])
    dmr<-unique(tmp$QTL)
    tmp<-subset(tmp, is.na(methStatus)!=1)
    tmp<-subset(tmp, expressionLog2>1)
    for(d in dmr){
      tmp2<-subset(tmp, QTL==d)
      if(length(unique(tmp2$methStatus))>=2){
        p1<- ggplot(tmp2, aes(x=as.factor(methStatus), y=expressionLog2, 
                              color=as.factor(methStatus))) +
          geom_boxplot(aes(color=as.factor(methStatus))) +
          scale_color_manual(values=c("#c33131", "#209f68")) +
          labs(title=genes[g], x="Methylation status", y ="Expression") +
          facet_grid(QTL ~ ., scales="free") +
          stat_compare_means(aes(label = paste0('p-value = ',..p.format..,' ',..p.signif..)),
                             cex=3,
                             method = "wilcox.test")
        p3<- ggplot(tmp2, aes(x=as.factor(methStatus), y=metaLevel, 
                              color=as.factor(methStatus))) +
          geom_boxplot(aes(color=as.factor(methStatus))) +
          scale_color_manual(values=c("#c33131", "#209f68")) +
          labs(title=genes[g], x="Methylation status", y ="Metabolite levels") +
          facet_grid(QTL ~ ., scales="free") +
          stat_compare_means(aes(label = paste0('p-value = ',..p.format..,' ',..p.signif..)),
                             cex=3, method = "wilcox.test")
        combined_plot <- grid.arrange(p1, p3, ncol = 1)
        ggsave(paste0(outPlot,'03.1_Meta-Expression-Meth_',genes[g],'_',d, "_correlation.png"), 
               combined_plot, width = 40, height = 20, units = "cm")
        ggsave(paste0(outPlot,'03.2_Meta-Expression-Meth_',genes[g],'_',d, "_correlation.pdf"), 
               combined_plot, width = 40, height = 20, units = "cm")
    }
  }
}

  ### get correlations and pvalues
  toSave<-c()
  for(g in 1:length(genes)){
    tmp<-subset(merged_data, geneName==genes[g])
    mQTL<-unique(tmp$QTL)
    tmp<-subset(tmp,expressionLog2>1)
    tmp<-subset(tmp, is.na(methStatus)!=1)
    for(l in mQTL){
        tmp2<-subset(tmp, QTL== l)
          nAlleles<-tmp2 %>% 
            group_by(methStatus) %>%
            summarise(count =n())
        if(dim(nAlleles)[1]>=2 & min(nAlleles$count)>=2){
          corPearson<-cor.test(tmp2$expressionLog2,tmp2$methLevels)
          pValExpMeth_Pearson<-corPearson$p.value
          corExpMeth_Pearson<-corPearson$estimate
      
          corSpearman<-cor.test(tmp2$expressionLog2,tmp2$methLevels, method = 'spearman')
          pValExpMeth_Spearman<-corSpearman$p.value
          corExpMeth_Spearman<-corSpearman$estimate
      
    } else {
      pValExpMeth_Pearson<-NA
      corExpMeth_Pearson<-NA
      pValExpMeth_Spearman<-NA
      corExpMeth_Spearman<-NA
    }
    tmp3<-cbind.data.frame(l,genes[g],corExpMeth_Pearson,pValExpMeth_Pearson,
                           corExpMeth_Spearman,pValExpMeth_Spearman)
    toSave<-rbind.data.frame(toSave, tmp3)
  }
}
  colnames(toSave)<-c('QTL','genes', 'corExpMeth_Pearson','pValExpMeth_Pearson',
                    'corExpMeth_Spearman','pValExpMeth_Spearman')
  head(toSave)
  sigCor<-toSave %>%
    filter(pValExpMeth_Pearson<=0.05)

  head(sigCor)
  merged_data[1:6,1:6]

  sigCor <- sigCor %>%
    left_join(merged_data, by = c("QTL" = "QTL", "genes" = "geneName")) %>%
    select("QTL", "kinship", "resultType","genes", "corExpMeth_Pearson", "pValExpMeth_Pearson","corExpMeth_Spearman", "pValExpMeth_Spearman",
           "metabolite","distance", "geneFunction") %>%
    distinct()

data.table::fwrite(sigCor, 
                   file= paste0(outDir,'06.1_sig-correlation_expression-meth.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")

allCor <- toSave %>%
  left_join(merged_data, by = c("QTL" = "QTL", "genes" = "geneName")) %>%
  select("QTL", "kinship", "resultType","genes", "corExpMeth_Pearson", "pValExpMeth_Pearson","corExpMeth_Spearman", "pValExpMeth_Spearman",
         "metabolite","distance", "geneFunction") %>%
  distinct()

data.table::fwrite(allCor, 
                   file= paste0(outDir,'06.2_all-correlation_expression-meth.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")
