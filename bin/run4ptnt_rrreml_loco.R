#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rrBLUP)
  library(data.table)
  library(dplyr)
  library(parallel)
  library(scales)
})

# -----------------------
# DIRS
# -----------------------
inDir      <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/ba_formated-data/"
clusterDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/bb_sample-clustering/"
outDir     <- "/mnt/disk2/vibanez/10_data-analysis/Fig5/ac_genomic_prediction/bc_predicted-values/"

# Crear subdirectorios necesarios
dir.create(file.path(outDir, "effects_ptnt"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outDir, "params_ptnt"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outDir, "summary_ptnt"), recursive = TRUE, showWarnings = FALSE)

# -----------------------
# ARGS
# -----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  mm      <- args[1]
  tt      <- args[2]
  tissue  <- args[3]
  dataSet <- args[4]
} else {
  # Defaults para testing local
  mm      <- "SNP"
  tt      <- "Tm003"
  tissue  <- "leaf"
  dataSet <- "snp-wo-na-all-dmrs"
}

# -----------------------
# LOAD DATA
# -----------------------
# Markers
markerP <- data.table::fread(
  paste0(inDir, mm, "_", tissue, "-metabolite.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE, fill = TRUE, na.strings = "NA"
)
sampleNames <- markerP$Samples
markerP$Samples <- NULL
if ("V1" %in% colnames(markerP)) markerP$V1 <- NULL
markerP <- as.data.frame(lapply(markerP, as.numeric))
rownames(markerP) <- sampleNames

# Traits
traits <- data.table::fread(
  paste0(inDir, "Traits_", tissue, "-metabolite.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE, fill = TRUE, na.strings = "NA"
)
rownames(traits) <- traits$V1
traits$V1 <- NULL
if (!(tt %in% colnames(traits))) {
  stop("Trait not found in traits table: ", tt)
}

# Construir y (log2) y alinear con X después
y_vec <- as.vector(traits[, tt])
names(y_vec) <- rownames(traits)

# Clusters (LOCO groups)
cluster_file <- paste0(clusterDir, "genetic-distance-clustering_cluster.tsv")
cluster_info <- data.table::fread(
  cluster_file, sep = "\t", header = TRUE, data.table = FALSE, fill = TRUE, na.strings = "NA"
)
colnames(cluster_info) <- c("Sample", "Cluster")

# -----------------------
# MATRICES X, y y alineación con foldid
# -----------------------
X <- as.matrix(markerP)
# Alinear y a X y transformar (log2)
y <- y_vec[rownames(X)]
if (any(is.na(y))) stop("NA in y after aligning to X; check sample overlap.")
y <- log2(y)

# Alinear foldid a X
foldid <- cluster_info$Cluster[ match(rownames(X), cluster_info$Sample) ]
if (any(is.na(foldid))) stop("Some samples in X not found in cluster file; check names.")
if (anyDuplicated(colnames(X)) > 0) stop("Duplicated marker/DMR column names in X.")

stopifnot(identical(rownames(X), names(y)))

set.seed(123)

# -----------------------
# HELPERS
# -----------------------
compute_marker_stats <- function(u, X_train, top_frac = 0.01) {
  stopifnot(length(u) == ncol(X_train))
  varZ <- apply(X_train, 2, var, na.rm = TRUE)
  VE   <- (u^2) * varZ
  propVE <- if (sum(VE, na.rm = TRUE) > 0) VE / sum(VE, na.rm = TRUE) else VE
  
  ord <- order(abs(u), decreasing = TRUE)
  rank_vec <- integer(length(u)); rank_vec[ord] <- seq_along(u)
  k <- max(1, round(top_frac * length(u)))
  top_idx <- ord[seq_len(k)]
  
  data.frame(
    Marker    = colnames(X_train),
    Effect    = as.numeric(u),
    AbsEffect = abs(u),
    Rank      = rank_vec,
    VE        = VE,
    PropVE    = propVE,
    TopK      = seq_len(length(u)) %in% top_idx,
    stringsAsFactors = FALSE
  )
}

perform_loco_cv_rrreml <- function(X, y, foldid, outDir, mm, tt, tissue, dataSet, ncores = 6) {
  cluster_ids <- unique(foldid)
  
  run_one <- function(cluster) {
    train_idx <- which(foldid != cluster)
    val_idx   <- which(foldid == cluster)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx, , drop = FALSE]
    y_val   <- y[val_idx]
    
    tryCatch({
      fit <- rrBLUP::mixed.solve(y_train, Z = X_train)
      
      y_pred <- as.numeric(X_val %*% fit$u + as.numeric(fit$beta[1]))
      PA   <- suppressWarnings(as.numeric(cor(y_val, y_pred, use = "complete.obs")))
      MSPE <- mean((y_val - y_pred)^2, na.rm = TRUE)
      MAPE <- mean(abs(y_val - y_pred), na.rm = TRUE)
      
      stats <- compute_marker_stats(fit$u, X_train, top_frac = 0.01)
      stats$Cluster <- cluster
      
      # Guardar efectos por fold
      dir.create(file.path(outDir, "effects_ptnt"), recursive = TRUE, showWarnings = FALSE)
      outfile <- paste0(outDir, "effects_ptnt/Effects_", mm, "_", tt, "_", tissue, "_", dataSet, "_Cluster", cluster, ".tsv")
      data.table::fwrite(stats[, c("Marker","Effect","AbsEffect","Rank","VE","PropVE","TopK")],
                         file = outfile, sep = "\t", quote = FALSE)
      
      # Guardar parámetros del modelo por fold
      beta0 <- as.numeric(fit$beta[1])
      Vu    <- as.numeric(fit$Vu)
      Ve    <- as.numeric(fit$Ve)
      params <- data.frame(
        Cluster = cluster, Intercept = beta0, Vu = Vu, Ve = Ve,
        nTrain = length(y_train), nVal = length(y_val),
        Trait = tt, Marker = mm, Tissue = tissue, Dataset = dataSet,
        stringsAsFactors = FALSE
      )
      dir.create(file.path(outDir, "params_ptnt"), recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(
        params,
        file = paste0(outDir, "params_ptnt/Params_", mm, "_", tt, "_", tissue, "_", dataSet, "_Cluster", cluster, ".tsv"),
        sep = "\t", quote = FALSE
      )
      
      list(
        perf  = data.frame(Cluster = cluster, PA = PA, MSPE = MSPE, MAPE = MAPE, Success = TRUE),
        stats = stats
      )
    }, error = function(e) {
      message("Error in cluster ", cluster, ": ", e$message)
      list(
        perf  = data.frame(Cluster = cluster, PA = NA, MSPE = NA, MAPE = NA, Success = FALSE),
        stats = NULL
      )
    })
  }
  
  res <- mclapply(cluster_ids, run_one, mc.cores = min(length(cluster_ids), ncores))
  perf_df  <- bind_rows(lapply(res, `[[`, "perf"))
  stats_df <- bind_rows(lapply(res, `[[`, "stats"))
  perf_df$Success_Rate <- mean(perf_df$Success, na.rm = TRUE)
  
  list(perf = perf_df, stats = stats_df)
}

# -----------------------
# RUN RR-REML LOCO
# -----------------------
res <- perform_loco_cv_rrreml(X, y, foldid, outDir, mm, tt, tissue, dataSet, ncores = 6)

# Save performance (as in your original naming)
perf <- res$perf
perf$Model        <- "RR-REML"
perf$Marker       <- mm
perf$Trait        <- tt
perf$Dataset_Type <- dataSet

data.table::fwrite(
  perf,
  file = paste0(outDir, "LOCO_", mm, "_", tt, "_", tissue, "_", dataSet, "_ptnt-version.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------
# CROSS-FOLD SUMMARY OF MARKERS
# -----------------------
all_stats <- res$stats
if (!is.null(all_stats) && nrow(all_stats) > 0) {
  summary_markers <- all_stats %>%
    group_by(Marker) %>%
    summarise(
      nFolds       = n(),
      freq_top1pct = mean(TopK, na.rm = TRUE),       # estabilidad
      mean_abs_u   = mean(AbsEffect, na.rm = TRUE),  # magnitud media
      mean_u       = mean(Effect, na.rm = TRUE),
      mean_rank    = mean(Rank, na.rm = TRUE),       # menor es mejor
      mean_propVE  = mean(PropVE, na.rm = TRUE),
      med_propVE   = median(PropVE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      score = rescale(freq_top1pct, to = c(0,1)) +
        0.7 * rescale(mean_abs_u,  to = c(0,1)) +
        0.7 * rescale(mean_propVE, to = c(0,1))
    ) %>%
    arrange(desc(score), mean_rank)
  
  data.table::fwrite(
    summary_markers,
    file = paste0(outDir, "summary_ptnt/MarkerSummary_", mm, "_", tt, "_", tissue, "_", dataSet, "_ptnt-version.tsv"),
    sep = "\t", quote = FALSE
  )
} else {
  warning("No marker stats to summarise (did all folds fail?).")
}

# -----------------------
# FINAL REFIT ON ALL DATA (modelo único para despliegue)
# -----------------------
message("Refit final con todas las muestras...")

# Guardar orden de marcadores (para reproducir en producción)
final_dir <- file.path(outDir, "final_model")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
marker_order <- colnames(X)
data.table::fwrite(
  data.frame(Marker = marker_order),
  file = file.path(final_dir, paste0("MarkersOrder_", mm, "_", tt, "_", tissue, "_", dataSet, "_ptnt-version.tsv")),
  sep = "\t", quote = FALSE
)

# Ajuste rrBLUP con todo el dataset
fit_full <- rrBLUP::mixed.solve(y, Z = X)
u_full   <- as.numeric(fit_full$u)        # pesos por DMR/marcador (w*)
mu_full  <- as.numeric(fit_full$beta[1])  # intercepto (mu*)

# 1) Guardar pesos por marcador (w*)
effects_full <- data.frame(
  Marker = colnames(X),
  Effect = u_full,
  stringsAsFactors = FALSE
)
data.table::fwrite(
  effects_full,
  file = file.path(final_dir, paste0("EffectsFull_", mm, "_", tt, "_", tissue, "_", dataSet, "_ptnt-version.tsv")),
  sep = "\t", quote = FALSE
)

# 2) Guardar intercepto y varianzas
params_full <- data.frame(
  Intercept = mu_full,
  Vu = as.numeric(fit_full$Vu),
  Ve = as.numeric(fit_full$Ve),
  n = nrow(X),
  p = ncol(X),
  Marker = mm, Trait = tt, Tissue = tissue, Dataset = dataSet,
  Timestamp = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
data.table::fwrite(
  params_full,
  file = file.path(final_dir, paste0("ParamsFull_", mm, "_", tt, "_", tissue, "_", dataSet, "_ptnt-version.tsv")),
  sep = "\t", quote = FALSE
)

# 3) (Opcional) Predicciones in-sample para control
y_hat_full <- as.numeric(X %*% u_full + mu_full)
pred_full <- data.frame(
  Sample = rownames(X),
  y_obs = y,
  y_hat = y_hat_full,
  stringsAsFactors = FALSE
)
data.table::fwrite(
  pred_full,
  file = file.path(final_dir, paste0("TrainPredFull_", mm, "_", tt, "_", tissue, "_", dataSet, "_ptnt-version.tsv")),
  sep = "\t", quote = FALSE
)

message("Listo: modelo final guardado. Ecuación: y_hat = mu_full + x^T u_full")
