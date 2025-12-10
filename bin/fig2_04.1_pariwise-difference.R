suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(ggpointdensity)
  library(viridis)
  library(segmented)
})
setwd("/mnt/disk2/vibanez/10_data-analysis/Fig2/ae_pairwise-divergences")
outDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/results/"
outPlot<-"/mnt/disk2/vibanez/10_data-analysis/Fig2/plots/"

treeDir<-"/mnt/disk2/vibanez/10_data-analysis/Fig1/"

treeOrder<- data.table::fread(paste0(treeDir,'general.SNPs.mdist.gz'), 
                                sep = '\t', data.table = F, 
                                fill = TRUE, na.string=c("NA"), nThread = 20) 
head(treeOrder)
colNames<-data.table::fread(paste0(treeDir, "general.SNPs.mdist.id"), 
                              sep = '\t', header = FALSE,
                              data.table = F, fill = TRUE, 
                              na.string=c("NA"), nThread = 20)[,2] 
  
  # rename samples to match
colNames<-gsub('^0_','',colNames)
colNames<-gsub('S_','S-',colNames)
colNames<-gsub('_1','',colNames)
colNames<-gsub('BGV006775','S-pimBGV006775',colNames)
colnames(treeOrder)<-colNames
rownames(treeOrder)<-colNames
head(treeOrder)
  
sampleTab<-data.table::fread(paste0(treeDir,"00_sampleTab.tsv"), 
                               sep = '\t', data.table = FALSE, fill = TRUE, 
                               na.string = "NA", nThread = 20)
sampleTab<-subset(sampleTab, Organ =="leaf")
head(sampleTab)
toKeep<-c('IID1','IID2','OBS_CT','DIFF_CT')
  # Load files
  
input<-list.files(pattern = ".sdiff.summary", full.names = FALSE)
processed_files <- list()
for ( i in input ){
    anno<-unlist(strsplit(i, '.', fixed=T))[1]
    mm<-unlist(strsplit(i, '.', fixed=T))[2]
    tmp <- data.table::fread(paste0(anno,'.',mm,'.',"sdiff.summary"), sep = '\t', data.table = T, fill = TRUE, 
                             na.string = "NA", nThread = 20)
    tmp<-tmp[,..toKeep]
    if(mm == "SNPs"){
      tmp$IID1<-gsub('_1','',gsub('S_','S-',tmp$IID1))
      tmp$IID2<-gsub('_1','',gsub('S_','S-',tmp$IID2))
      tmp<-tmp[!(tmp$IID1 %in% 'S-penLA0716') & !(tmp$IID2 %in% 'S-penLA0716'),]
      tmp$IID1<-gsub('BGV006775','S-pimBGV006775',tmp$IID1)
      tmp$IID2<-gsub('BGV006775','S-pimBGV006775',tmp$IID2)
      tmp<-tmp[order(tmp$IID1, tmp$IID2),]
      tmp$comp <- ifelse(tmp$IID1 < tmp$IID2, paste0(tmp$IID1, ':', tmp$IID2), paste0(tmp$IID2, ':', tmp$IID1))
      tmp$DIFF_CT<-tmp$DIFF_CT/1000000 #1e6
    } else {
      tmp$IID1<-gsub('_leaf','',tmp$IID1)
      tmp$IID2<-gsub('_leaf','',tmp$IID2)
      tmp$IID1<-gsub('LA2838A','S-lycLA2838A',tmp$IID1)
      tmp$IID2<-gsub('LA2838A','S-lycLA2838A',tmp$IID2)
      tmp<-tmp[!(tmp$IID1 %in% "S-lycLA2845") & !(tmp$IID2 %in% "S-lycLA2845"),]
      tmp<-tmp[order(tmp$IID1, tmp$IID2),]
      tmp$comp <- ifelse(tmp$IID1 < tmp$IID2, paste0(tmp$IID1, ':', tmp$IID2), paste0(tmp$IID2, ':', tmp$IID1))
      if(mm=="CG-DMR"){
        scale<-100000
      }else{
        scale<-1000000
      }
      tmp$DIFF_CT<-tmp$DIFF_CT/scale #1e4
      dim(tmp)
    }
    setnames(tmp, old = c("OBS_CT", "DIFF_CT"), new = c(paste0("n", mm),mm))
    processed_files[[mm]] <- tmp
}
merged_data <- Reduce(function(x, y) merge(x, y, by = c("comp"), all = TRUE), processed_files)
merged_data<-merged_data[,c(1,2,3,12,13,4,5,8,9)]
colnames(merged_data)<-c('sampleComparison','acc1','acc2','nSNP','SNPs',
                           'nC-DMR','C-DMR','nCG-DMR','CG-DMR')
head(merged_data)
merged_data<-subset(merged_data, acc1!='S-pimLA1578' & acc2!='S-pimLA1578')
merged_data<-subset(merged_data, acc1!='TS-244' & acc2!='TS-244')

## add group information
group<-unique(sampleTab$Group)
pairwiseDifferences<-c()
for( g in group){
    tmpTab<-subset(sampleTab, Group==g)
    tmpDF<-merged_data[merged_data$acc1 %in% tmpTab$BiseqName,]
    tmpDF$group1<-g
    group2 <- c()  # Initialize an empty vector for group2
    for (a in 1:nrow(tmpDF)) {
      a2<-as.character(tmpDF[a,3])
      g2 <- subset(sampleTab, BiseqName == a2)$Group
      if(length(g2)>1){
        print(a2)
      }
      group2 <- c(group2, g2)
    }
    tmpDF$group2 <- group2
    pairwiseDifferences<-rbind.data.frame(pairwiseDifferences, tmpDF)
}

pairwiseDifferences$groupComparison<-paste0(pairwiseDifferences$group1,":",pairwiseDifferences$group2)
  
## add SNPs distance from the phylogenetic three
pairwiseDifferences$SNPdistance <- apply(pairwiseDifferences, 1, function(x) {
    treeOrder[x['acc1'], x['acc2']]})
  
  keep<-c("sampleComparison","acc1", "acc2","SNPdistance","SNPs","C-DMR", "CG-DMR","group1", 
          "group2", "groupComparison")
pairwiseDifferences<-pairwiseDifferences[,..keep]

#########################################
### get plots
# get color with all samples
color_map <- c('PIM:PIM' = "#00a681ff", 'SLC:SLC' = "#ed83b5ff", 'SLL:SLL' = "#00a4deff",
               'SLC:PIM' = "#eaed83", 'PIM:SLC' = "#eaed83", 'SLC:SLL' = "#a900de",
               'SLL:SLC' = "#a900de", 'SLL:PIM' = "#de3a00", 'PIM:SLL' = "#de3a00",
               'WILD:WILD' = "grey")
rmIntrogressed<-c("S-lycEA00371","S-lycEA00940","S-lycEA00990","S-lycEA01155",
                  "S-lycEA02054","S-lycLA1090","S-lycLA2838A", "S-lycLYC1410",
                  "S-lycLYC1969", "S-lycLYC3153", "TS-10","TS-112","TS-122","TS-128",
                  "TS-138", "TS-140","TS-212","TS-3","TS-518","TS-530","TS-540", "TS-554",
                  "TS-560", "TS-566","TS-578", "TS-579","TS-637","TS-638","TS-656","TS-67",
                  "TS-8","TS-9","TS-95")

correlation<-cor(pairwiseDifferences$SNPdistance, pairwiseDifferences$SNPs)

pairwiseDifferences <- pairwiseDifferences %>%
  mutate(highlight = ifelse(acc1 %in% rmIntrogressed | acc2 %in% rmIntrogressed, "Introgression", "Normal"))

ggplot(pairwiseDifferences, aes(x = SNPdistance, y = SNPs, color = groupComparison, shape = highlight)) +
  geom_point(alpha = 0.5, aes(size = highlight)) +  # Scatter points with different sizes for highlight
  geom_smooth(aes(x = SNPdistance, y = SNPs), method = "lm", color = "blue", se = TRUE, inherit.aes = FALSE) +  # Single regression line
  labs(
    title = paste("Correlation between SNPdistance and # SNPs: ", round(cor(pairwiseDifferences$SNPdistance, pairwiseDifferences$SNPs, use = "complete.obs"), 2)),
    x = "SNP Distance",
    y = "Pairwise SNP Difference"
  ) +
  scale_color_manual(values = color_map) +  # Apply color map to groupComparison
  scale_shape_manual(values = c("Introgression" = 17, "Normal" = 16)) +  # Different shapes for highlight and normal
  scale_size_manual(values = c("Introgression" = 3, "Normal" = 1)) +  # Larger points for highlighted samples
  theme_minimal()

# Save the plot
ggsave(paste0(outPlot, "04.0_SNPdistance-numbSNPs_highlight.pdf"), width = 20, height = 20, units = "cm")

#
# Data labels for the plots
labels <- list(
  SNPdistance_CG_DMR = list(x ="Genetic distance (1-IBS)", 
                            y = expression(paste("pairwise CG-DMRs differences (", 10^5, ")"))),
  SNPdistance_C_DMR = list(x = "Genetic distance (1-IBS)", 
                           y = expression(paste("pairwise C-DMRs differences (", 10^6, ")")))
)

# Function to create and save the plot
create_boxplot_with_dots <- function(data, x_var, y_var, color_by, output_file, x_label, y_label) {
  breaks <- seq(0, 1, by = 0.025)  # Fixed range from 0 to 1 with 0.025 intervals
  data$bin <- cut(data[[x_var]], breaks = breaks, include.lowest = TRUE, right = TRUE)
  # Format bin labels
  bin_labels <- levels(data$bin)
  formatted_labels <- sub("\\((.+),(.+)\\]", "\\1-\\2", bin_labels)  # Format like '0-0.025'
  data$bin <- factor(data$bin, levels = bin_labels, labels = formatted_labels)
  
  # Select bins to display on x-axis
  display_bins <- formatted_labels[seq(1, length(formatted_labels), by = 2)]  # Show every 4th bin
  display_bins <- c(display_bins, formatted_labels[length(formatted_labels)])  # Ensure last bin is shown
  
  # Plot
  p <- ggplot(data, aes(x = bin, y = .data[[y_var]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot without color
    geom_jitter(aes(color = groupComparison), width = 0.1, alpha = 0.4, size = 0.5) +  # Color only the dots
    xlab(x_label) +
    ylab(y_label) +
    scale_color_manual(values = color_map) +  # Apply color map only to dots
    scale_x_discrete(breaks = display_bins) +  # Show selected bins on x-axis
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
  
  # Save the plot
  ggsave(paste0(outPlot, output_file), plot = p, width = 20, height = 20, units = "cm")
}
# Example usage
create_boxplot_with_dots(pairwiseDifferences, "SNPdistance", "CG-DMR", 
                         "groupComparison", 
                         "04.1.CG-DMR_SNP_distance_allSamples_boxplot.pdf", 
                         labels$SNPdistance_CG_DMR$x, labels$SNPdistance_CG_DMR$y)
create_boxplot_with_dots(pairwiseDifferences, "SNPdistance", "C-DMR", 
                         "groupComparison", 
                         "04.1.C-DMR_SNP_distance_allSamples_boxplot.pdf", 
                         labels$SNPdistance_C_DMR$x, labels$SNPdistance_C_DMR$y)
  
#### remove samples
woWild<-subset(pairwiseDifferences, group1!='WILD' & group2!='WILD')
head(woWild)
##### remove introgressions
# samples with more than 1Mb of introgressions:  subset(summ, maxIntro>=1000000)$acc 

woWildwoIntro<-woWild[!(woWild$acc1 %in% rmIntrogressed),]
woWildwoIntro<-woWildwoIntro[!(woWildwoIntro$acc2 %in% rmIntrogressed),]
woWildwoIntro

#https://itecnotes.com/tecnote/r-piecewiseregression-with-r-plotting-the-segments/
  
# Function to create segmented plots
create_segmented_plot <- function(data, x_var, y_var, color_by, output_file, x_label, y_label) {
    x_var_sym <- sym(x_var)
    y_var_sym <- sym(y_var)
    
    # Perform segmented regression
    seg_model <- segmented(
      lm(as.formula(paste0("`", y_var, "` ~ 0 + `", x_var, "`")), data = data),
      seg.Z = as.formula(paste0("~ `", x_var, "`"))
    )
    
    predicted <- data.frame(x = data[[x_var]], y = broken.line(seg_model)$fit)
    colnames(predicted) <- c(x_var, y_var)
    predicted <- predicted[order(predicted[[x_var]], predicted[[y_var]]), ]
    
    # Calculate differences and derivatives
    predicted$diff_x <- c(0, diff(predicted[[x_var]]))
    predicted$diff_y <- c(0, diff(predicted[[y_var]]))
    predicted$derivative <- predicted$diff_y / predicted$diff_x
    abrupt_changes <- which(abs(diff(predicted$derivative * 100)) > 1)
    
    if(length(abrupt_changes) == 0){
      xChange <- NA
      xRate <- NA
    } else {
      abrupt_changes <- min(abrupt_changes)
      xChange <- predicted[[x_var]][abrupt_changes]
      xRate <- predicted$derivative[abrupt_changes]
    }
    
    # Handle missing y values for annotation
    max_y <- max(data[[y_var]], na.rm = TRUE)
    if (is.infinite(max_y)) {
      max_y <- 0
    }
    
    p <- ggplot(data, aes(x = !!x_var_sym, y = !!y_var_sym)) +
        geom_point(aes(color = groupComparison), alpha = 0.4, size = 5) +
        geom_line(data = predicted, aes(x = !!x_var_sym, y = !!y_var_sym), color = 'blue') +
        xlab(x_label) +
        ylab(y_label) +
        scale_color_manual(values = color_map) +
        theme_bw()
      
    if(!is.na(xChange)){
        p <- p + geom_vline(xintercept = xChange, linetype = "dashed", color = "red") +
          annotate("text", x = xChange, y = max_y * 0.8, 
                   label = paste("rho:", round(xRate, 2)), 
                   vjust = 1, hjust = 0.5, color = "red", size = 3)
    ggsave(paste0(outPlot,output_file), plot = p, width = 20, height = 20, units = "cm")
    }
  }
create_segmented_plot(woWildwoIntro, "SNPdistance", "CG-DMR", "groupComparison",
                        "04.2.CG-DMR_SNP_distance_clean_colcomparison.pdf", 
                        labels$SNPdistance_CG_DMR$x, labels$SNPdistance_CG_DMR$y)
  
create_segmented_plot(woWildwoIntro, "SNPdistance", "C-DMR", "groupComparison",
                        "04.2.C-DMR_SNP_distance_clean_colcomparison.pdf", 
                        labels$SNPdistance_C_DMR$x, labels$SNPdistance_C_DMR$y)
###
create_combined_plot <- function(data, x_var, y_var, color_by, output_file, x_label, y_label) {
  x_var_sym <- sym(x_var)
  y_var_sym <- sym(y_var)
  
  # Define bins for the x-axis
  breaks <- seq(0, 1, by = 0.025)  # Fixed range from 0 to 1 with 0.025 intervals
  data$bin <- cut(data[[x_var]], breaks = breaks, include.lowest = TRUE, right = TRUE)
  
  # Format bin labels
  bin_labels <- levels(data$bin)
  formatted_labels <- sub("\\((.+),(.+)\\]", "\\1-\\2", bin_labels)
  data$bin <- factor(data$bin, levels = bin_labels, labels = formatted_labels)
  
  # Aggregate data by bins for fitting
  aggregated <- data %>%
    group_by(bin) %>%
    summarize(
      mean_x = mean(.data[[x_var]], na.rm = TRUE),
      mean_y = mean(.data[[y_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Fit loess model to aggregated data
  loess_fit <- loess(mean_y ~ mean_x, data = aggregated)
  aggregated$loess_y <- predict(loess_fit, aggregated$mean_x)
  
  # Find the maximum point (where the line stops increasing)
  derivatives <- diff(aggregated$loess_y) / diff(aggregated$mean_x)
  max_point <- which(derivatives <= 0)[1]  # First point where the rate of increase stops
  
  if (!is.na(max_point) && !is.null(max_point)) {
    xMax <- aggregated$mean_x[max_point]
    yMax <- aggregated$loess_y[max_point]
  } else {
    xMax <- NA
    yMax <- NA
  }
  
  # Handle missing y values
  max_y <- max(data[[y_var]], na.rm = TRUE)
  if (is.infinite(max_y)) {
    max_y <- 0
  }
  
  # Plot
  p <- ggplot(data, aes(x = bin, y = .data[[y_var]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = groupComparison), width = 0.1, alpha = 0.4, size = 0.5) +  # Color only the dots
    geom_line(data = aggregated, aes(x = bin, y = loess_y), color = 'blue', group = 1, inherit.aes = FALSE) +
    xlab(x_label) +
    ylab(y_label) +
    scale_color_manual(values = color_map) +  # Apply color map only to dots
    scale_x_discrete(breaks = formatted_labels[seq(1, length(formatted_labels), by = 4)]) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (!is.na(xMax)) {
    p <- p +
      geom_vline(xintercept = xMax, linetype = "dashed", color = "red") +
      annotate("text", x = xMax, y = max_y * 0.9, 
               label = paste("Max Point:", round(xMax, 3)), 
               vjust = -0.5, hjust = 1, color = "red", size = 3)
  }
  
  # Ensure output directory exists
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Save the plot
  tryCatch({
    ggsave(paste0(outPlot, output_file), plot = p, width = 20, height = 20, units = "cm")
    print(paste("Plot saved successfully at", output_file))
  }, error = function(e) {
    print(paste("Error during ggsave:", e$message))
  })
}

create_combined_plot(woWildwoIntro, "SNPdistance", "CG-DMR", "groupComparison", 
                     "04.3.CG-DMR_SNP_distance_combined_plot.pdf", 
                     labels$SNPdistance_CG_DMR$x, labels$SNPdistance_CG_DMR$y)
create_combined_plot(woWildwoIntro, "SNPdistance", "C-DMR", "groupComparison", 
                     "04.3.C-DMR_SNP_distance_combined_plot.pdf", 
                     labels$SNPdistance_C_DMR$x, labels$SNPdistance_C_DMR$y)

create_combined_plot(pairwiseDifferences, "SNPdistance", "CG-DMR", "groupComparison", 
                     "04.4.CG-DMR_SNP_allSamples_distance_combined_plot.pdf", 
                     labels$SNPdistance_CG_DMR$x, labels$SNPdistance_CG_DMR$y)
create_combined_plot(pairwiseDifferences, "SNPdistance", "C-DMR", "groupComparison", 
                     "04.4.C-DMR_SNP_allSamples_distance_combined_plot.pdf", 
                     labels$SNPdistance_C_DMR$x, labels$SNPdistance_C_DMR$y)
###
#remove only introgressed lin

woIntro<-pairwiseDifferences[!(pairwiseDifferences$acc1 %in% rmIntrogressed),]
woIntro<-pairwiseDifferences[!(pairwiseDifferences$acc2 %in% rmIntrogressed),]
woIntro
create_combined_plot(woIntro, "SNPdistance", "CG-DMR", "groupComparison", 
                     "04.5.CG-DMR_SNP_woIntro_distance_combined_plot.pdf", 
                     labels$SNPdistance_CG_DMR$x, labels$SNPdistance_CG_DMR$y)
create_combined_plot(woIntro, "SNPdistance", "C-DMR", "groupComparison", 
                     "04.5.C-DMR_SNP_woIntro_distance_combined_plot.pdf", 
                     labels$SNPdistance_C_DMR$x, labels$SNPdistance_C_DMR$y)


