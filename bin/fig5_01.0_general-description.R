{
suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(eulerr)
  library(ggridges)
  library(ggplot2)
})
outDir<-"10_data-analysis/Fig5/ab_data-analysis/results/"
outPlot<-"10_data-analysis/Fig5/ab_data-analysis/plots/"
#### get numbers of metabolites associ-non-associated
basedir<-"10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/"
asso<-c('sig/','nonSig/','nonSig-GIF/')
out<-c()
betaVals<-c()
QTLs<-c()
for(a in asso){
  resultType<-gsub('/','',a)
  input<-list.files(path=paste0(basedir,a), pattern = ".QTL", full.names = FALSE)
  for(i in input){
    metabolite<-unlist(strsplit(i, '.', fixed=TRUE))[1]
    mm<-unlist(strsplit(i, '_', fixed=TRUE))[2]
    kinship<-unlist(strsplit(i, '_', fixed=TRUE))[3]
    
    index<-unlist(strsplit(i, '_', fixed=TRUE))[4]
    index<-gsub('.QTL','',index)
    
    h2<- data.table::fread(gsub('.QTL','.reml',paste0(basedir,a,i)), 
                           sep = '\t', data.table = FALSE, fill = TRUE, 
                           na.string=c("NA"), nThread = 20)
    H2<-h2[4,1]/(h2[4,1]+h2[5,1])
    h2<-h2[6,1]
    
    if(a =='nonSig/'){
      nLoci<-0
      mmEffect<-0
    } else {
        df<- data.table::fread(paste0(basedir,a,i), sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("NA"), nThread = 20)
        nLoci<-nrow(df)
        mmEffect<-df$V2
        qtl<-df$V1
      }
      tmp<-cbind.data.frame(metabolite,mm,kinship,index,resultType,nLoci, h2, H2)
      tmpBeta<-cbind.data.frame(metabolite,mm,kinship,index,resultType, mmEffect)
      tmpQTL<-cbind.data.frame(metabolite,mm,kinship,index,resultType, qtl)
      out<-rbind.data.frame(out, tmp)
      betaVals<-rbind.data.frame(tmpBeta, betaVals)
      QTLs<-rbind.data.frame(QTLs, tmpQTL)
    }
}
data.table::fwrite(out, file= paste0(outDir,'01.0_general-information.tsv'), quote=F,
                       row.names=F,col.names = T,sep="\t")
data.table::fwrite(QTLs, file= paste0(outDir,'01.1_QTLs.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")
data.table::fwrite(betaVals, file= paste0(outDir,'01.2_beta-values.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")

out <- out[order(out$metabolite), ]

head(out)
# plot beta values
{
  betaVals$absValues<-abs(betaVals$mmEffect)
  betaVals$kinship <- factor(betaVals$kinship, levels = c("SNP", "DMR", "DMR-SNP"))
  bSplit<-betaVals %>%
    filter(resultType=="sig" & index=='aBN' & mm!="DMR-SNP")
  ggplot(bSplit, aes(x=mm, y=absValues, fill=mm, color=mm)) +
    geom_boxplot(alpha=0.1, outlier.shape = 16, outlier.size = 1.5) +  # Remove redundant color and set outliers
    geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) + 
    scale_color_manual(values=c("#c06a80", "#999999"))+
    scale_fill_manual(values=c("#c06a80","#999999"))+
    labs(x = "Marker",y = "Absolute Beta Values") +
    facet_wrap(.~kinship)+
    theme(
      plot.background = element_rect(fill = "white", color = "black"),  # Add border
      #panel.grid.major = element_line(color = "gray85"),  # Adjust gridline color
      axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(outPlot,"01.1_beta_absValues.pdf"), width = 30, height = 20, units = "cm")
  
  ggplot(bSplit, aes(x=mm, y=mmEffect, fill=mm, color=mm)) +
    geom_boxplot(alpha=0.1, outlier.shape = 16, outlier.size = 1.5) + 
    geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +  # Add jittered point
    scale_color_manual(values=c("#c06a80", "#999999"))+
    scale_fill_manual(values=c("#c06a80","#999999"))+
    facet_wrap(.~kinship)+
    labs(x = "Marker",y = "Beta Values") +
    theme(plot.background = element_rect(fill = "white", color = "black"),  # Add border
          #panel.grid.major = element_line(color = "gray85"),  # Adjust gridline color
          axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
          legend.position = "none")
  ggsave(paste0(outPlot,"01.2_beta_values.pdf"), width = 30, height = 20, units = "cm")
  }
# plot H2
{

  H2split<-out %>%
  filter(mm!="DMR-SNP" & index=="aBN" & kinship!="DMR-SNP")

  summary_df<-H2split %>%
    group_by(mm ,kinship) %>%
    dplyr::summarise(mean_h2 =mean(h2, na.rm=T))

    # Create the plot with specified facet order
  ggplot(H2split, aes(x = h2, color=kinship, fill=kinship)) + 
    geom_density(alpha=0.1) +
    labs(title="", x="Heritability", y = "Density") +
    scale_color_manual(values=c("#c06a80", "#999999"))+
    scale_fill_manual(values=c("#c06a80","#999999"))+
    theme_ridges(font_size = 10, center_axis_labels = T) +
    theme_minimal() +
    geom_vline(data = summary_df, aes(xintercept = mean_h2, color = mm), 
               linetype = "dashed", linewidth = 0.8) +
    theme(plot.background = element_rect(fill = "white"))
  ggsave(paste0(outPlot,"02.1_Heritability_individual-mm.pdf"), width = 30, height = 20, units = "cm")
  
  h2<-out %>%
    filter(mm!="DMR-SNP"  & kinship=="DMR-SNP")
  
  summary_df<-h2 %>%
    group_by(mm, resultType) %>%
    dplyr::summarise(mean_h2 =mean(h2, na.rm=T))
  
  # Create the plot with specified facet order
  ggplot(h2, aes(x = h2, color = mm, fill = mm)) + 
    geom_density(alpha = 0.1) +
    labs(title = "", x = "Heritability", y = "Density") +
    scale_color_manual(values = c("#c06a80", "#999999")) +
    scale_fill_manual(values = c("#c06a80", "#999999")) +
    theme_ridges(font_size = 10, center_axis_labels = TRUE) +
    facet_grid(. ~ resultType, 
               labeller = labeller(
                 resultType = c(
                   "sig" = "Significant", 
                   "nonSig-GIF" = "Sig with GI", 
                   "NonSig" = "Non-Significant"
                 )
               )) +
    theme_minimal() +
    geom_vline(
      data = summary_df, 
      aes(xintercept = mean_h2, color = mm), 
      linetype = "dashed", linewidth = 0.8
    ) +
    theme(plot.background = element_rect(fill = "white"))
  ggsave(paste0(outPlot,"02.2_Heritability_DMR-SNP.pdf"), width = 30, height = 20, units = "cm")
}
# number of meta with QTL
summary<-out %>%
  filter(mm!="DMR-SNP" & index =="aBN") %>%
  group_by(mm, kinship, resultType) %>%
  dplyr::summarise(nMetabolites =n())

data.table::fwrite(summary, 
                   file= paste0(outDir,'01.3_number-metabolite-associated.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")
#### get numbers of loci
summLoci<-out %>%
  group_by(mm,kinship, resultType) %>%
  dplyr::summarize(nLoci=sum(nLoci))
data.table::fwrite(summLoci, 
                   file= paste0(outDir,'01.4_number-loci-metabolite_associated.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")


