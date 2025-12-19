suppressPackageStartupMessages({
  library(glmnet)
  library(rrBLUP)
  library(data.table)
  library(dplyr)
  library(rbridge)
  library(ncvreg)
  library(sparsenet)
  library(parallel)
  library(dplyr)
})

# Define directories
inDir <- "/10_data-analysis/Fig5/ac_genomic_prediction/ba_formated-data/"
clusterDir<- "10_data-analysis/Fig5/ac_genomic_prediction/bb_sample-clustering/"
outDir <- "10_data-analysis/Fig5/ac_genomic_prediction/bc_predicted-values/"

# Arguments
args <- commandArgs(trailingOnly = TRUE)
mm <- args[1]  # Marker type: "SNP" or "DMR"
tt <- args[2]  # Trait
tissue <- args[3]  # Tissue
dataSet <- args[4]  # Cluster type: "SNP", "DMR", "Combined"

# Load marker data
markerP <- data.table::fread(paste0(inDir, mm, '_', tissue, '-metabolite.tsv'), 
                             sep = '\t', header = TRUE, data.table = FALSE, fill = TRUE, na.strings = "NA")
sampleNames <- markerP$Samples
markerP$Samples <- NULL
markerP <- as.data.frame(lapply(markerP, as.numeric))
rownames(markerP) <- sampleNames

# Load trait data
traits <- data.table::fread(paste0(inDir,'Traits_', tissue, '-metabolite.tsv'), 
                            sep = '\t', header = TRUE, data.table = FALSE, fill = TRUE, na.strings = "NA")
rownames(traits) <- traits$V1
traits$V1 <- NULL
y <- log2(as.vector(traits[, tt]))

# Load cluster information based on the specified cluster type
#cluster_file <- paste0(clusterDir, mm, "_cluster_",dataSet,".tsv") # 4 supplementary data
cluster_file <- paste0(clusterDir, "genetic-distance-clustering_cluster.tsv")
cluster_info <- data.table::fread(cluster_file, sep = '\t', header = TRUE, data.table = FALSE, fill = TRUE, na.strings = "NA")
colnames(cluster_info) <- c("Sample", "Cluster")
cluster_info <- cluster_info[match(rownames(markerP), cluster_info$Sample), ]
foldid <- cluster_info$Cluster

# Function for LOCO CV
perform_loco_cv <- function(model_type, X, y, foldid) {
  cluster_ids <- unique(foldid)
  run_cluster <- function(cluster) {
    train_indices <- which(foldid != cluster)
    val_indices <- which(foldid == cluster)
    
    X_train <- X[train_indices, ]
    y_train <- y[train_indices]
    X_val <- X[val_indices, ]
    y_val <- y[val_indices]
    
    tryCatch({
      if (model_type == "RR-REML") {
        model <- mixed.solve(y_train, Z = X_train)
        y_pred <- X_val %*% model$u + model$beta[1]
      } else if (model_type == "RR-CV") {
        model <- cv.glmnet(X_train, y_train, family = "gaussian", type = "mse", alpha = 0, nlambda = 100, seed = 123)
        y_pred <- predict(model, X_val, s = "lambda.min")
      } else if (model_type == "LASSO") {
        model <- cv.ncvreg(X_train, y_train, penalty = "lasso", seed = 123)
        y_pred <- predict(model, X_val, s = "lambda.min")
      } else if (model_type == "ENET") {
        model <- cv.glmnet(X_train, y_train, family = "gaussian", alpha = 0.5, nlambda = 100, seed = 123)
        y_pred <- predict(model, X_val, s = "lambda.min")
      } else if (model_type == "Ridge") {
        model <- cv.glmnet(X_train, y_train, family = "gaussian", alpha = 0, nlambda = 100, seed = 123)
        y_pred <- predict(model, X_val, s = "lambda.min")
      } else if (model_type == "sENET") {
        model <- cv.sparsenet(X_train, y_train, trace.it = FALSE)
        y_pred <- predict(model, X_val, which = "parms.min")
      }
      
      # Metrics
      PA <- as.numeric(cor(y_val, y_pred, use = "complete.obs"))
      MSPE <- as.numeric(mean((y_val - y_pred)^2, na.rm = TRUE))
      MAPE <- as.numeric(mean(abs(y_val - y_pred), na.rm = TRUE))
      
      if (is.na(PA)) {
        return(data.frame(Cluster = cluster, PA = NA, MSPE = NA, MAPE = NA, Success = FALSE))
      } else {
        return(data.frame(Cluster = cluster, PA = PA, MSPE = MSPE, MAPE = MAPE, Success = TRUE))
      }
      
    }, error = function(e) {
      cat("Error in cluster", cluster, "for model", model_type, ":", e$message, "\n")
      return(data.frame(Cluster = cluster, PA = NA, MSPE = NA, MAPE = NA, Success = FALSE))
    })
  }
  
  results <- mclapply(cluster_ids, run_cluster, mc.cores = min(length(cluster_ids), 10))
  performance <- do.call(rbind, results)
  
  # Add success rate
  success_rate <- sum(performance$Success) / length(cluster_ids)
  performance$Success_Rate <- success_rate
  
  # Optional: drop Success column
  performance$Success <- NULL
  
  return(performance)
}
# Models to evaluate
models <- c("RR-REML","RR-CV", "Ridge", "LASSO","ENET","sENET")
# Run LOCO CV for each model
n_cores <- 6
all_results <- mclapply(models, function(model) {
  cat("Running model:", model, "\n")
  results <- perform_loco_cv(model, as.matrix(markerP), y, foldid)
  results$Model <- model
  results$Marker <- mm
  results$Trait <- tt
  results$Dataset_Type <- dataSet
  return(results)
}, mc.cores = n_cores)

names(all_results) <- models

# Combine results into a single data frame
final_results <- do.call(rbind, all_results)

write.table(final_results, file = paste0(outDir, "LOCO_", mm, "_", tt, "_", tissue, "_", dataSet, ".tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
