suppressPackageStartupMessages({
  library(data.table)
})
wd<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/"

# merge h2
check<-c('sig/','nonSig/','nonSig-GIF')
out <- list()
for (c in check) {
  input <- list.files(path = file.path(wd, c), pattern = ".reml", full.names = TRUE)
  result <- rbindlist(lapply(input, function(file_path) {
    df <- fread(file_path, sep = '\t', fill = TRUE, na.string = "NA", nThread = 20)
    DMR <- unlist(strsplit(basename(file_path), ".", fixed = TRUE))[1]
    type<-gsub('/','',c)
    df <- transpose(df)
    df[, DMR := DMR]
    df[, type := type]
    #print(df)
    return(df)
  }))
  out[[c]] <- result
}
h2merged <- rbindlist(out, use.names = TRUE, fill = TRUE)
colnames(h2merged)<-c("logVariance","log","delta","sigmaGen","sigmaError", "h2",'DMR','type')
data.table::fwrite(h2merged, file= paste0(outDir,'00.1_merged-H2_allDMRs.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")
# merge beta values
input <- list.files(path = file.path(wd, 'sig'), pattern = ".mQTL", full.names = TRUE)
betaOut <- rbindlist(lapply(input, function(file_path) {
    df <- fread(file_path, sep = '\t', fill = TRUE, na.string = "NA", nThread = 20)
    DMR <- unlist(strsplit(basename(file_path), ".", fixed = TRUE))[1]
    df[, DMR := DMR]
    return(df)
}))
betaOut
colnames(betaOut)<-c("SNPid","beta","betaSD","pvalues","DMRid")
data.table::fwrite(betaOut, file= paste0(outDir,'00.1_merged-beta-values_sig-DMRs.tsv'), quote=F,
                   row.names=F,col.names = T, sep="\t")
