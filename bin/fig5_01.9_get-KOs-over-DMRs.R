suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

inDir<-"/mnt/disk2/vibanez/09_KO-processing/natural-experimental_DMR/"
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig5/ab_data-analysis/results/"
# data from Camus:/mnt/disk2/vibanez/natural-experimental_DMR
#################### LOAD DATA TO MERGE  ######################################################
input<-list.files(path=inDir,pattern = "bed", full.names = FALSE)

replace_ones <- function(x) {ifelse(x == "1", 0, x)}
replace_one_with_text <- function(x) {gsub("^1,.*", "1", x)}

out <- c()
for (i in input){
  dmr<-unlist(strsplit(i, '.', fixed=TRUE))[1]
  dmr<-gsub('ko-','',dmr)

  dt<-data.table::fread(paste0(inDir,i), sep = '\t',
                        data.table = T, fill = TRUE, check.names=FALSE,
                        na.string=c("NA"), nThread = 10)
  dt$V10<-NULL
  dt <- dt %>%
    mutate(across(-c(chr, start, end), replace_ones)) %>%
    mutate(across(-c(chr, start, end), replace_one_with_text))
  dt$DMRtarget<-dmr
  out<-rbind.data.frame(out, dt)
}
head(out)
out <- out %>%
  data.table::melt(id.vars = c("chr", "start", "end", "DMRtarget"), 
                 variable.name = "ko", value.name = "target")%>%
  filter(target==1)
out$DMR<-paste0(gsub('SL2.50','',out$chr),"_",out$DMRtarget,'_',out$start)

collapsed_data <- out %>%
  group_by(chr, start, end,DMR) %>%
  summarize(
    ko = paste(unique(ko), collapse = ","),
    DMRtarget = paste(unique(DMRtarget), collapse = ",")
  ) %>%
  ungroup()

data.table::fwrite(collapsed_data, 
                   file= paste0(outDir,'07.1_DMRs-over-KO.tsv'), quote=F,
                   row.names=F,col.names = T,sep="\t")
