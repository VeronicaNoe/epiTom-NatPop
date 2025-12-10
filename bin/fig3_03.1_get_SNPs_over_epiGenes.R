#!/home/vibanez/anaconda3/envs/mKit/bin/Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(qqman)
  library(ggplot2)
  library(ggrepel)
})


outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/plots/"

input<-list.files(path=outDir,pattern = '.closest-gene', full.names = FALSE)

# Read and process files
dt <- lapply(input, function(input) {
  dmr <- unlist(strsplit(basename(input), "_", fixed = TRUE))[2]
  dmr <- unlist(strsplit(dmr, ".", fixed = TRUE))[1]
  out <- fread(paste0(outDir,input), sep = '\t', header = F, fill = TRUE, na.strings = "")
  out$dmr <- dmr
  return(out)
})
dt <- rbindlist(dt)

colnames(dt)<-c('chrSNP', 'startSNP','endSNP','chisq','pvalue','nPvalues',
                'chrGene','startGene','endGene','strand','toRemove',
                'geneName','geneFunction','toRemove','DMR')

toSort<-c('chrSNP', 'startSNP','endSNP','DMR','chisq','pvalue','nPvalues',
          'chrGene','startGene','endGene','strand','geneName','geneFunction')

dt<-dt[,..toSort]
length(unique(dt$geneName))
data.table::fwrite(dt, file= paste0(outDir,'03.0_SNPs_over_genes.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")

#dt_filter <- dt[pvalue !=0]
# check some gens
#dt_filter[grepl("Solyc10g005130", geneName)]

######3 load list of genes in epi pathways
geneList <- fread("epigenetic_path_genes.tsv", sep = '\t', header = F, fill = TRUE, na.strings = "")
out<-c()
for( g in geneList$V1){
  tmp<-dt[grepl(g, geneName)]
  if(dim(tmp)[1]>0){
    gList<-geneList %>%
      filter(V1==g)
    tmp$Name<-paste0(gList$V2, collapse = ',')
  }
  out<-rbind.data.frame(out, tmp)
}
out[, snpList := paste0(sprintf("%02d", chrSNP), ":", startSNP)]
out[, .(
  DMR = paste0(unique(DMR), collapse = ","),
  nSNPsOverEpiGenes = .N,
  minPvalue = min(pvalue),
  unique_geneFunction = paste(unique(geneFunction), collapse = "; "),
  unique_name = paste(unique(Name), collapse = "; ")
), by = geneName]


snpList<-unique(out$snpList)
length(snpList)
data.table::fwrite(as.data.table(snpList), file= paste0(outDir,'03.1_SNPs_over_epiGenes_list.tsv'), quote=F,
                   row.names=F,col.names = F, sep="\t")

data.table::fwrite(out, file= paste0(outDir,'03.2_SNPs_over_epiGenes.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")
# do Camus the ~/bin/check_top5chisq.sh
# Group by 'DMR' and summarize
summary_dt <- out[, .(
  DMR = paste0(unique(DMR), collapse = ","),
  nSNPsOverEpiGenes = .N,
  minPvalue = min(pvalue),
  unique_geneFunction = paste(unique(geneFunction), collapse = "; "),
  unique_name = paste(unique(Name), collapse = "; ")
), by = geneName]

data.table::fwrite(summary_dt, file= paste0(outDir,'03.3_SNPs_over_epiGenes_summary.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")

#nrpde=IPR006592 Solyc01g096390
#ddm1 = IPR000330 Solyc02g085390.2
#cmt3=IPR001525 Solyc01g006100 5?
#cmt3=IPR001525 Solyc12g100330
#met1= IPR017198 Solyc11g030600 1copy
#kyp=IPR003105  Solyc02g094520 9?
