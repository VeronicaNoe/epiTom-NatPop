suppressPackageStartupMessages({
  require("data.table")
#  library(plyr)
library(dplyr)
})
inDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/be_check-GIF/cc_lambda-tables"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/"
###########################
## A GIF range of 0.975 to 1.025 allows for a maximum deviation of 2.5% from
## the null distribution of test statistics in either direction.
## This means that your analysis can tolerate a maximum of 2.5% inflation or deflation
## of test statistics from what would be expected under the null hypothesis.
lambda <- data.table(dmrPos = character(), dmr = character(), nAssociation = integer(),
                         lambda = character(), sigLambda = character())
input<-list.files(path=inDir, pattern = ".tsv", full.names = FALSE)
t<-0
for (i in input){
  print(length(input) - t)
  dmrPos<-unlist(strsplit(i, ".", fixed = TRUE))[1]
  dmrPos<-gsub("_lambda", "",dmrPos)
  df<-data.table::fread(paste0(inDir,"/",i), sep = '\t',data.table = T, fill = T,
                         header = T, na.string=c("NA"), nThread = 20)
  df[, dmrPos := dmrPos]
  df[, sigLambda := ifelse(lambda >= 0.975 & lambda <= 1.025, "linked",'linked_w_GI')]
  lambda <- rbindlist(list(lambda, df), fill = TRUE)
  t <- t +1
}
data.table::fwrite(lambda, file=paste0(outDir,"00.0_results-adjuted-lambda.tsv"), quote=F,
                   row.names=F,col.names = T,sep="\t")

nonSig<-subset(lambda, sigLambda=="linked_w_GI")$dmrPos
data.table::fwrite(as.data.table(nonSig), file=paste0(outDir,"00.0_list2move-from-sig2nonsig.tsv"),
                   quote=F, row.names=F,col.names = F,sep="\t")
lambda %>%
  group_by(dmr, sigLambda) %>%
  dplyr::summarize(count = n())
sig<-subset(lambda, sigLambda=="linked")$dmrPos
data.table::fwrite(as.data.table(sig), file=paste0(outDir,"00.0_list_of_sigAssociations.tsv"),
                   quote=F, row.names=F,col.names = T,sep="\t")
