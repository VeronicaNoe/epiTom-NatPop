suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
})
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig2/ab_get-LD")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/plots/"

# zeroDist<-c()
# for(i in 1:length(input)){
#   df<- data.table::fread(input[i], sep = '\t', data.table = T, 
#                          fill = TRUE, na.string=c("NA"), nThread = 20, check.names = FALSE) 
#   df<-subset(df, distance<=0)
#   print(input[i])
#   print(head(df))
#   if(dim(df)[1]==0){next}
#   #info to add
#   sample<-unlist(strsplit(input[i],".", fixed = TRUE))[1]
#   #create each vector
#   ANNOTATION<-unlist(strsplit(sample,"_", fixed = TRUE))[1]
#   #<-rep(unlist(strsplit(sample,"_", fixed = TRUE))[3],times=nrow(df))
#   if(is.na(unlist(strsplit(sample,"_", fixed = TRUE))[2])==TRUE){
#     DMR<-rep("SNPs",times=nrow(df))
#   } else{
#     DMR<-rep(unlist(strsplit(sample,"_", fixed = TRUE))[2],times=nrow(df))
#   }
#   tmp<-cbind.data.frame(df,  ANNOTATION, GROUP)
#   colnames(tmp)[3]<-"avg_R2"
#   zeroDist<-rbind.data.frame(zeroDist, tmp)
#   colnames(zeroDist)[3]<-"avg_R2"
# }
# head(zeroDist)

input<-list.files(pattern = ".ld_decay", full.names = FALSE)
input<-grep('bins', input, invert = T, value = T)
ld<-c()
for(i in 1:length(input)){
  df <- data.table::fread(input[i], sep = '\t', data.table = TRUE, 
                          fill = TRUE, na.string=c("NA"), 
                          nThread = 20, check.names = FALSE)  # Adjust encoding if necessary
  df[, chr := sub("^b'", "", chr)]
  df[, chr := sub("'$", "", chr)]
  #info to add
  marker<-unlist(strsplit(input[i],".", fixed = TRUE))[1]
  marker<-unlist(strsplit(marker,"_", fixed = TRUE))[3]
  ann<-unlist(strsplit(input[i],"_", fixed = TRUE))[2]
  group<-unlist(strsplit(input[i],"_", fixed = TRUE))[1]
  df$group<-group
  df$annotation<-ann
  df$marker<-marker
  #tmp<-cbind.data.frame(df, DMR, CHR, ANNOTATION,GROUP)
  ld<-rbind.data.frame(ld, df)
}
sort(unique(ld$chr))
min(ld$distance)

ld$distanceKb<-ld$distance/1000 #kbases
# average by chromosomes
summ<-plyr::ddply(ld, c("distanceKb", "marker", "group"), summarise,
                        mean = mean(R2, na.rm=TRUE),
                        sd= mean(std, na.rm=TRUE),
                        sem = sd(R2,na.rm=TRUE)/sqrt(length(R2)))
tmp<-subset(summ, group=="allAcc")
ggplot(tmp, aes(x=distanceBin, y=mean, group=group)) +
  geom_line(aes(color=marker), linewidth=1)+
  scale_x_continuous(limits = c(0,14), breaks = seq(0,14, by=1))+
  xlab("Distance (Kb)") + ylab(expression(italic(r)^2))#+
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, group=DMR), 
  #             alpha = 0.1, color = "grey") +
  #scale_color_manual(values=c("#820a86", "#ffca7b", "#999999"))+
  #facet_grid(~group)
ggsave(paste0(outDir,"02.1.LD-decay.png"),  width = 20, height = 20, units = "cm")

# only general 
toPlot<-subset(ld, GROUP=="general")
summ<-plyr::ddply(toPlot, c("distanceBin","GROUP", "DMR"), summarise,
                  mean = mean(avg_R2, na.rm=TRUE),
                  sd= mean(std, na.rm=TRUE),
                  sem = sd(avg_R2,na.rm=TRUE)/sqrt(length(avg_R2)))
summ[summ==0.5]<-0.0
# plot LD decay
ggplot(summ, aes(x=distanceBin, y=mean, group=DMR)) +
  geom_line(aes(color=DMR), linewidth=1)+
  scale_x_continuous(limits = c(0,15), breaks = seq(0,30, by=1))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by=0.2))+
  xlab("Distance (Kb)") + ylab(expression(italic(r)^2))+
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, group=DMR), 
  #             alpha = 0.1, color = "grey") +
  scale_color_manual(values=c("#820a86", "#ffca7b", "#999999"))+
  facet_grid(~GROUP)
ggsave(paste0(outDir,"02.2.all_general_LD-decay.png"),  width = 20, height = 20, units = "cm")

# split by chromosomes
# chr<-unique(ld$chr)
# for(c in chr){
#   toPlot<-subset(ld, chr==c & GROUP=="all")
#   # plot LD decay
#   ggplot(toPlot, aes(x=distance, y=avg_R2, group=DMR)) +
#     geom_line(aes(color=DMR), linewidth=1)+
#     xlim(0, 12)+
#     xlab("Distance (Kb)") + ylab(expression(italic(r)^2))+
#     scale_color_manual(values=c("#820a86", "#ffca7b", "#999999"))+
#     facet_wrap(~GROUP)
#   ggsave(paste0(outDir,"02.3.",c,"_LD-decay.png"),  width = 20, height = 20, units = "cm")
# }
