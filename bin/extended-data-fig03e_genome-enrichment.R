#!/usr/bin/env Rscript

# Extended Data Fig. 3e: genome-wide enrichment by distance to the closest TE.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ============================================================================
# 1. Paths and settings -- edit here when running locally
# ============================================================================
outDir  <- "/mnt/disk2/vibanez/zenodo/ExtDataFig03/"
outPlot <- "/mnt/disk2/vibanez/zenodo/ExtDataFig03/"
allDataFile   <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/results/02.07_epiallele_sample.tsv"
closestTEFile <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/results/00_allGenes_closestTE.tsv"

prefix <- "gene-epiallele_genome-enrichment"
classes <- c("gbM", "teM", "UM")
classColours <- c(gbM = "#ffca7b", teM = "#820a86", UM = "grey50")
expectedColour <- "grey30"
binWidth <- 500L
maxDistance <- 3000L
nRandom <- 1000L
randomSeed <- 20261017L

# Hide Outside only in the figures, to match the supplied panel.
# Outside is ALWAYS retained in the universe, random draws and test family.
showOutsideBin <- FALSE
# "stars" keeps the supplied presentation. "adjustedP" prints numeric adjusted
# P values above bins instead. Full numeric/scientific results are always saved.
significanceStyle <- "stars"
plotWidthCm <- 20
plotHeightCm <- 8.5
yDivisor <- 1000
barOffset <- 0.20
barWidth <- 0.36

stopifnot(binWidth > 0L, maxDistance > 0L,
          maxDistance %% binWidth == 0L,
          nRandom >= 10L, nRandom == as.integer(nRandom),
          significanceStyle %in% c("stars", "adjustedP"))
if (getRversion() < "3.6.0") stop("Use R >= 3.6.0 for reproducible sampling.")
if (packageVersion("ggplot2") < "3.4.0") stop("Use ggplot2 >= 3.4.0.")
for (f in c(allDataFile, closestTEFile)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}
for (d in unique(c(outDir, outPlot))) {
  if (!dir.exists(d) && !dir.create(d, recursive = TRUE)) {
    stop("Cannot create output directory: ", d)
  }
}

# ============================================================================
# 2. Helpers
# ============================================================================
require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}
write_tsv <- function(x, suffix) {
  fwrite(x, file.path(outDir, paste0(prefix, suffix)), sep = "\t",
         quote = FALSE, na = "NA", col.names = TRUE)
}
format_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# Maintain log probabilities for very small P values: an underflowed numeric
# zero is not reported as an exactly zero P value. Use the same probability-
# ordering definition and relative tolerance as R's two-sided fisher.test.
fisher_log_p <- function(observed, binTotal, universeTotal, classTotal) {
  support <- seq.int(max(0L, classTotal - (universeTotal - binTotal)),
                     min(classTotal, binTotal))
  lp <- dhyper(support, binTotal, universeTotal - binTotal, classTotal, log = TRUE)
  lp_observed <- dhyper(observed, binTotal, universeTotal - binTotal,
                       classTotal, log = TRUE)
  selected <- lp[lp <= lp_observed + log1p(1e-7)]
  peak <- max(selected)
  min(0, peak + log(sum(exp(selected - peak))))
}
bh_log <- function(lp) {
  # BH adjusted probabilities in log space, including underflowed P values.
  m <- length(lp)
  o <- order(lp, decreasing = TRUE)
  adjusted <- pmin(0, cummin(lp[o] + log(m) - log(m:1L)))
  answer <- numeric(m)
  answer[o] <- adjusted
  answer
}
format_log_p <- function(lp, digits = 12L) {
  vapply(lp, function(z) {
    if (is.na(z)) return("NA")
    if (z == 0) return("1")
    if (!is.finite(z)) return("numerical_underflow")
    exponent <- floor(z / log(10))
    mantissa <- exp(z - exponent * log(10))
    if (round(mantissa, digits) >= 10) {
      mantissa <- mantissa / 10
      exponent <- exponent + 1
    }
    sprintf(paste0("%.", digits, "fe%+d"), mantissa, as.integer(exponent))
  }, character(1))
}
star_label <- function(logp) {
  fifelse(logp < log(0.001), "***",
          fifelse(logp < log(0.01), "**", fifelse(logp < log(0.05), "*", "")))
}

# ============================================================================
# 3. Gene universe and mutually exclusive class assignments
# ============================================================================
allData <- fread(allDataFile, sep = "\t", header = TRUE,
                 select = c("geneName", "geneEpiallele"),
                 na.strings = c("NA", ""), check.names = FALSE)
allGeneTE <- fread(closestTEFile, sep = "\t", header = TRUE,
                   na.strings = c("NA", ""), check.names = FALSE)
require_columns(allData, c("geneName", "geneEpiallele"), "Epiallele table")
require_columns(allGeneTE,
                c("geneName", "geneStrand", "teChr", "teStart", "teEnd", "signedDistance"),
                "Closest-TE table")
allData[, `:=`(geneName = as.character(geneName),
               geneEpiallele = as.character(geneEpiallele))]
allGeneTE[, geneName := as.character(geneName)]
if (anyNA(allGeneTE$geneName) || any(allGeneTE$geneName == "")) {
  stop("Every row of the gene universe must have a geneName.")
}
if (anyDuplicated(allGeneTE$geneName)) {
  stop("The closest-TE table must contain exactly one row per gene. Resolve ties upstream.")
}
if (!nrow(allGeneTE)) stop("The closest-TE gene universe is empty.")

classMap <- unique(allData[geneEpiallele %chin% classes & !is.na(geneName),
                           .(geneName, geneEpiallele)])
conflicts <- classMap[, .(nClasses = uniqueN(geneEpiallele)), by = geneName][nClasses > 1L]
if (nrow(conflicts)) {
  stop("Multiple class assignments for genes: ",
       paste(head(conflicts$geneName, 10L), collapse = ", "))
}
missingGenes <- setdiff(classMap$geneName, allGeneTE$geneName)
if (length(missingGenes)) {
  stop("Classified genes absent from the closest-TE universe (do not silently drop them): ",
       paste(head(missingGenes, 10L), collapse = ", "))
}
# Sort once so the same input rows in a different order give the same draws.
geneData <- copy(allGeneTE)
geneData[, geneClass := "Other"]
geneData[classMap, on = .(geneName), geneClass := i.geneEpiallele]
setorder(geneData, geneName)
for (column in c("teStart", "teEnd", "signedDistance")) {
  original <- geneData[[column]]
  numericValue <- suppressWarnings(as.numeric(original))
  invalid <- !is.na(original) & original != "." & is.na(numericValue)
  if (any(invalid)) stop("Non-numeric entries in ", column, "; fix the closest-TE table.")
  set(geneData, j = column, value = numericValue)
}

# ============================================================================
# 4. Preserve the original signed-distance bins and Outside category
# ============================================================================
up <- seq.int(-maxDistance, -binWidth, by = binWidth)
down <- seq.int(binWidth, maxDistance, by = binWidth)
allKeys <- c(as.character(up), "Overlap", as.character(down), "Outside")
labels <- setNames(c(sprintf("%.1f", up / 1000), "Overlap",
                     sprintf("%.1f", down / 1000), ">3 kb / no TE"), allKeys)
binDefinitions <- data.table(
  distanceKey = allKeys,
  distanceLabel = unname(labels),
  lowerBp = c(up, 0, down - binWidth, NA_real_),
  upperBp = c(up + binWidth, 0, down, NA_real_),
  interval = c(sprintf("[%d,%d)", up, up + binWidth), "{0}",
               sprintf("(%d,%d]", down - binWidth, down),
               "abs(signedDistance) > 3000 OR no TE hit")
)
geneData[, noTEHit := is.na(teChr) | teChr == "." | is.na(teStart) | teStart < 0]
if (geneData[!noTEHit & !is.finite(signedDistance), .N]) {
  stop("A gene with a TE hit has missing/non-finite signedDistance.")
}
geneData[, distanceKey := fcase(
  noTEHit, "Outside",
  abs(signedDistance) > maxDistance, "Outside",
  signedDistance == 0, "Overlap",
  default = as.character(sign(signedDistance) * binWidth *
                         ceiling(abs(signedDistance) / binWidth))
)]
if (anyNA(geneData$distanceKey) || any(!geneData$distanceKey %chin% allKeys)) {
  stop("Unrecognized distance category; check closest-TE distances.")
}
# Same active-bin rule as the original script: omit entirely empty universe bins.
activeKeys <- allKeys[allKeys %chin% geneData$distanceKey]
universeTotal <- nrow(geneData)
binTotals <- tabulate(match(geneData$distanceKey, activeKeys), nbins = length(activeKeys))
nClass <- vapply(classes, function(cl) sum(geneData$geneClass == cl), integer(1))
if (any(nClass == 0L) || any(nClass == universeTotal)) {
  stop("Each requested class must have >0 and <N genes in the universe. Counts: ",
       paste(classes, nClass, sep = "=", collapse = "; "))
}

write_tsv(geneData, "_gene-universe.tsv")
binDefinitions[, universeGenes := as.integer(table(factor(geneData$distanceKey, levels = allKeys)))]
write_tsv(binDefinitions, "_bin-definitions.tsv")
sampleSizes <- data.table(geneClass = classes, targetTotal = nClass,
                          referenceTotal = universeTotal - nClass,
                          universeTotal = universeTotal,
                          nRandomSets = nRandom, seed = randomSeed,
                          nBinsTested = length(activeKeys))
write_tsv(sampleSizes, "_sample-sizes.tsv")
print(sampleSizes)

# ============================================================================
# 5. Exact per-bin tests against all non-class genes; not against Other alone
# ============================================================================
testList <- vector("list", length(classes) * length(activeKeys))
index <- 0L
for (className in classes) {
  k <- nClass[[className]]
  observed <- tabulate(match(geneData[geneClass == className, distanceKey], activeKeys),
                       nbins = length(activeKeys))
  for (j in seq_along(activeKeys)) {
    K <- binTotals[j]
    x <- observed[j]
    # Rows: class / non-class. Columns: bin / all remaining bins.
    contingency <- matrix(c(x, k - x, K - x, universeTotal - k - K + x),
                          nrow = 2L, byrow = TRUE)
    if (K == universeTotal) {
      # Degenerate bin: all genes are here, so there is no contrast.
      oddsRatio <- lower95OR <- upper95OR <- NA_real_
      pValue <- 1
      logP <- 0
      testStatus <- "Degenerate: all universe genes in this bin"
    } else {
      ft <- fisher.test(contingency, alternative = "two.sided", conf.level = 0.95)
      oddsRatio <- unname(ft$estimate)
      lower95OR <- ft$conf.int[1L]
      upper95OR <- ft$conf.int[2L]
      pValue <- ft$p.value
      logP <- if (pValue > 0) log(pValue) else fisher_log_p(x, K, universeTotal, k)
      testStatus <- "OK"
    }
    index <- index + 1L
    testList[[index]] <- data.table(
      comparison = "Genome", reference = "All non-class genes",
      geneClass = className, distanceKey = activeKeys[j],
      distanceCategory = unname(labels[activeKeys[j]]),
      targetTotal = k, referenceTotal = universeTotal - k,
      universeTotal = universeTotal, binTotal = K,
      targetInBin = x, targetOutsideBin = k - x,
      referenceInBin = K - x,
      referenceOutsideBin = universeTotal - k - K + x,
      observedCount = x, expectedMean = k * K / universeTotal,
      expectedLower95 = qhyper(0.025, K, universeTotal - K, k),
      expectedUpper95 = qhyper(0.975, K, universeTotal - K, k),
      oddsRatio = oddsRatio, lower95OR = lower95OR, upper95OR = upper95OR,
      pValue = pValue, logPValue = logP, testStatus = testStatus
    )
  }
}
results <- rbindlist(testList)
results[, logAdjustedP := bh_log(logPValue), by = geneClass]
results[, `:=`(
  adjustedP = exp(logAdjustedP),
  pValueScientific = format_log_p(logPValue),
  adjustedPScientific = format_log_p(logAdjustedP),
  pValueUnderflow = pValue == 0,
  foldEnrichment = observedCount / expectedMean,
  significance = star_label(logAdjustedP),
  direction = fcase(oddsRatio > 1, "Enriched", oddsRatio < 1, "Depleted",
                    default = "No difference"),
  nRandomSets = nRandom,
  multiplicity = "BH across all active bins within each gene class"
)]
write_tsv(results, "_per-bin_tests.tsv")

# ============================================================================
# 6. Actual random gene sets, used ONLY to display the null distribution
# ============================================================================
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(randomSeed)
universeBins <- match(geneData$distanceKey, activeKeys)
randomList <- vector("list", length(classes))
for (cc in seq_along(classes)) {
  className <- classes[cc]
  k <- nClass[[className]]
  message("Sampling ", nRandom, " random sets of ", k, " genes for ", className)
  counts <- matrix(0L, nrow = nRandom, ncol = length(activeKeys))
  for (draw in seq_len(nRandom)) {
    chosen <- sample.int(universeTotal, size = k, replace = FALSE)
    counts[draw, ] <- tabulate(universeBins[chosen], nbins = length(activeKeys))
  }
  # One sampled set supplies counts for ALL bins; never draw each bin separately.
  stopifnot(all(rowSums(counts) == k),
            all(sweep(counts, 2L, binTotals, FUN = "<=") | counts == 0L))
  randomList[[cc]] <- data.table(
    geneClass = className, iteration = rep(seq_len(nRandom), times = length(activeKeys)),
    distanceKey = rep(activeKeys, each = nRandom),
    randomCount = as.integer(counts), setSize = k, universeTotal = universeTotal
  )
}
randomCounts <- rbindlist(randomList)
write_tsv(randomCounts, "_random-set_counts.tsv")

# Explicit box summaries match the unweighted geom_boxplot type-7 quartiles.
boxSummary <- randomCounts[, {
  quartiles <- as.numeric(quantile(randomCount, c(0.25, 0.5, 0.75), type = 7))
  IQR <- quartiles[3L] - quartiles[1L]
  whiskerLow <- min(randomCount[randomCount >= quartiles[1L] - 1.5 * IQR])
  whiskerHigh <- max(randomCount[randomCount <= quartiles[3L] + 1.5 * IQR])
  .(nRandomSets = .N, randomMean = mean(randomCount), randomSD = sd(randomCount),
    minimum = min(randomCount), q25 = quartiles[1L], median = quartiles[2L],
    q75 = quartiles[3L], maximum = max(randomCount),
    lowerWhisker = whiskerLow, upperWhisker = whiskerHigh,
    nOutliers = sum(randomCount < whiskerLow | randomCount > whiskerHigh))
}, by = .(geneClass, distanceKey)]
write_tsv(boxSummary, "_boxplot_summary.tsv")

# ============================================================================
# 7. Plot data: preserve all tests even when Outside is hidden visually
# ============================================================================
plotKeys <- if (showOutsideBin) activeKeys else setdiff(activeKeys, "Outside")
if (!length(plotKeys)) stop("No genes fall in the displayed bins. Set showOutsideBin=TRUE.")
facetLabels <- setNames(paste0(classes, "\nn = ", format_n(nClass), " genes"), classes)
plotResults <- copy(results[distanceKey %chin% plotKeys])
plotRandom <- copy(randomCounts[distanceKey %chin% plotKeys])
plotBoxes <- copy(boxSummary[distanceKey %chin% plotKeys])
for (objectName in c("plotResults", "plotRandom", "plotBoxes")) {
  object <- get(objectName)
  object[, `:=`(geneClass = factor(geneClass, levels = classes),
                xIndex = match(distanceKey, plotKeys))]
  assign(objectName, object)
}
plotRandom[, expectedX := xIndex + barOffset]
plotResults[, `:=`(observedX = xIndex - barOffset,
                   expectedX = xIndex + barOffset)]
randomMax <- plotRandom[, .(maxRandomCount = max(randomCount)), by = .(geneClass, distanceKey)]
plotResults <- merge(plotResults, randomMax, by = c("geneClass", "distanceKey"), all.x = TRUE)
plotResults[, geneClass := factor(geneClass, levels = classes)]
plotResults[, panelMaximum := max(c(observedCount, expectedUpper95, maxRandomCount)), by = geneClass]
plotResults[, starY := pmax(observedCount, expectedUpper95, maxRandomCount) +
                       0.055 * pmax(panelMaximum, 1)]
plotResults[, annotation := if (significanceStyle == "stars") significance else
              paste0("adj.P=", format_log_p(logAdjustedP, digits = 1L))]

#significanceStyle <- "adjustedP" , "stars"
caption <- paste0("Universe: ", format_n(universeTotal), " genes; ", format_n(nRandom),
                  " random gene sets/class. Random draws are not biological replicates.\n",
                  if (significanceStyle == "stars")
                    "BH-adjusted P: * <0.05; ** <0.01; *** <0.001. " else "",
                  "Exact test results: accompanying per-bin Source Data.")
common_layers <- function() {
  list(
    facet_wrap(~geneClass, nrow = 1L, scales = "free_y",
               labeller = as_labeller(facetLabels)),
    scale_x_continuous(breaks = seq_along(plotKeys),
                       labels = unname(labels[plotKeys]),
                       limits = c(0.5, length(plotKeys) + 0.5), expand = c(0, 0)),
    scale_y_continuous(labels = label_number(accuracy = 0.01),
                       expand = expansion(mult = c(0, 0.13))),
    labs(x = "Signed distance to closest TE (kb)",
         y = "Number of genes (x 1,000)", tag = "e", caption = caption),
    theme_classic(base_size = 8, base_family = "sans"),
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 8),
          legend.position = "top", legend.title = element_blank(),
          legend.key.size = grid::unit(0.30, "cm"),
          panel.spacing = grid::unit(0.35, "cm"),
          plot.caption = element_text(size = 6, hjust = 0),
          plot.subtitle = element_text(size = 7),
          plot.tag = element_text(face = "bold", size = 11),
          plot.margin = margin(5, 5, 5, 5))
  )
}

# ============================================================================
# 8. Version A: observed bars, exact expected mean + 95% NULL interval,
#    and all individual random-set counts overlaid on the expected bars
# ============================================================================
pBar <- ggplot() +
  geom_col(data = plotResults,
           aes(x = observedX, y = observedCount / yDivisor, fill = geneClass),
           width = barWidth) +
  geom_col(data = plotResults,
           aes(x = expectedX, y = expectedMean / yDivisor, fill = "Expected"),
           width = barWidth) +
  geom_point(data = plotRandom, aes(x = expectedX, y = randomCount / yDivisor),
             position = position_jitter(width = 0.10, height = 0, seed = randomSeed + 1L),
             colour = "grey70", alpha = 0.16, size = 0.35, show.legend = FALSE) +
  geom_errorbar(data = plotResults,
                aes(x = expectedX, ymin = expectedLower95 / yDivisor,
                    ymax = expectedUpper95 / yDivisor),
                width = 0.13, linewidth = 0.30, colour = "black") +
  geom_text(data = plotResults,
             aes(x = xIndex, y = starY / yDivisor, label = annotation),
             size = if (significanceStyle == "stars") 2.4 else 1.8) +
  scale_fill_manual(values = c(classColours, Expected = expectedColour),
                     breaks = c(classes, "Expected"), drop = FALSE) +
  common_layers() +
  labs(subtitle = "Colored bars: observed; grey bars: expected; dots: random sets; whiskers: 95% null interval")

# ============================================================================
# 9. Version B: boxplots of random-set COUNTS, observed count as one diamond
# ============================================================================
pBox <- ggplot() +
  geom_boxplot(data = plotBoxes,
                aes(x = xIndex, ymin = lowerWhisker / yDivisor,
                    lower = q25 / yDivisor, middle = median / yDivisor,
                    upper = q75 / yDivisor, ymax = upperWhisker / yDivisor,
                    group = interaction(geneClass, distanceKey)),
                stat = "identity", orientation = "x", width = 0.55,
                fill = "grey85", colour = "grey30", linewidth = 0.30) +
  # All draws, including outliers, are shown here; outliers are not drawn twice.
  geom_point(data = plotRandom, aes(x = xIndex, y = randomCount / yDivisor),
              position = position_jitter(width = 0.16, height = 0, seed = randomSeed + 2L),
              colour = "grey50", alpha = 0.12, size = 0.35) +
  geom_point(data = plotResults,
              aes(x = xIndex, y = observedCount / yDivisor, fill = geneClass),
              shape = 23, size = 1.9, stroke = 0.35, colour = "black",
              show.legend = FALSE) +
  geom_text(data = plotResults,
              aes(x = xIndex, y = starY / yDivisor, label = annotation),
              size = if (significanceStyle == "stars") 2.4 else 1.8) +
  scale_fill_manual(values = classColours, drop = FALSE) +
  common_layers() +
  labs(subtitle = "Grey boxes and dots: random gene sets; colored diamonds: observed gene counts")

pdfDevice <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
barFile <- file.path(outPlot, paste0(prefix, "_observed-expected_counts.pdf"))
boxFile <- file.path(outPlot, paste0(prefix, "_random-set_boxplots.pdf"))
ggsave(barFile, pBar, device = pdfDevice, width = plotWidthCm,
       height = plotHeightCm, units = "cm", bg = "white")
ggsave(boxFile, pBox, device = pdfDevice, width = plotWidthCm,
       height = plotHeightCm, units = "cm", bg = "white")

# ============================================================================
# 10. Legends with actual n, complete definitions and exact per-bin P values
# ============================================================================
nText <- paste(paste0(classes, ": n = ", format_n(nClass), " genes"), collapse = "; ")
baseLegend <- paste0(
  "e, Distribution of gbM, teM and UM genes by signed distance to the closest TE ",
  "in 500-bp bins; Overlap denotes distance zero. ", nText, ". ",
  "The sampling universe contained N = ", format_n(universeTotal),
  " genes from xx00_allGenes_closestTE.tsv. This table is assumed to represent ",
  "the complete annotated-gene universe. Distance signs follow this input table; ",
  "strand orientation is not recomputed. For each class, ", format_n(nRandom),
  " random gene sets of the same size were drawn from the full universe, without ",
  "replacement within each draw; all genes were eligible again in the next draw ",
  "(seed = ", randomSeed, "). Random draws are computational null samples, not ",
  "biological or technical replicates. The counted unit is a unique gene. ",
  "Genes >", maxDistance / 1000, " kb from their closest TE or lacking a TE hit ",
  "were included in the universe, remaining-bin counts and multiple-testing family",
  if (showOutsideBin) "; they are displayed as Outside. " else ", but their Outside bin is not displayed. "
)
statLegend <- paste0(
  "Enrichment/depletion was assessed using two-sided Fisher exact tests on ",
  "2 x 2 tables (class versus all non-class genes; bin versus all remaining bins). ",
  "Benjamini-Hochberg adjustment was applied across the ", length(activeKeys),
  " active distance bins separately within each class. These are per-bin tests, ",
  "not an omnibus test, and are not tests across random draws. ",
  "Asterisks indicate adjusted P < 0.05 (*), < 0.01 (**) and < 0.001 (***). ",
  "Exact test P values, adjusted P values, odds ratios with conditional 95% ",
  "confidence intervals and all four contingency-table counts are supplied in ",
  prefix, "_per-bin_tests.tsv. Numeric probabilities may underflow; the ",
  "pValueScientific/adjustedPScientific and log columns retain their magnitude."
)
barLegend <- paste0(baseLegend,
  "Colored bars show the observed aggregate gene counts (one count per class/bin). ",
  "Grey bars show the analytical null expectation k*K/N, where k is the class ",
  "size and K is the number of universe genes in the bin. Dots show all random-set ",
  "counts. Error bars show the exact 2.5th and 97.5th hypergeometric percentiles ",
  "for random-set counts. They are pointwise null intervals (at least 95% coverage ",
  "because counts are discrete), NOT SD, SEM or confidence intervals of the mean. ",
  statLegend)
boxLegend <- paste0(baseLegend,
  "Grey boxplots summarize the ", format_n(nRandom),
  " random-set counts per class/bin: centre line, median; box limits, 25th and ",
  "75th percentiles; whiskers, smallest/largest random-set counts within 1.5 times ",
  "the interquartile range below/above the box limits. All random-set counts, ",
  "including those beyond the whiskers, are shown as horizontally jittered dots ",
  "without vertical jitter. Colored diamonds show the single observed aggregate ",
  "count for each class/bin, not a separate sample distribution. ", statLegend)
exactLegendLines <- results[, paste0(geneClass, ", bin ", distanceCategory,
                                     ": P = ", pValueScientific,
                                     "; BH-adjusted P = ", adjustedPScientific)]
writeLines(c("BAR VERSION", barLegend, "", "BOXPLOT VERSION", boxLegend, "",
             "EXACT PER-BIN TEST VALUES (for a detailed legend or Source Data)",
             exactLegendLines),
           file.path(outDir, paste0(prefix, "_legends.txt")))

settings <- data.table(
  parameter = c("allDataFile", "closestTEFile", "universeTotal", "classes", "binWidth", "maxDistance",
                "nRandom", "randomSeed", "showOutsideBin", "significanceStyle", "BH_family", "universe_rule"),
  value = c(allDataFile, closestTEFile, universeTotal, paste(classes, collapse = ";"), binWidth, maxDistance,
            nRandom, randomSeed, showOutsideBin, significanceStyle,
            "All active bins, including Outside, separately within each class",
            "All rows in closest-TE input; exactly one row per gene; no distance restriction")
)
write_tsv(settings, "_settings.tsv")
writeLines(capture.output(sessionInfo()), file.path(outDir, paste0(prefix, "_sessionInfo.txt")))

message("Done.\nBar plot: ", barFile, "\nBoxplot: ", boxFile,
        "\nTests, random counts, gene universe, sample sizes and legends: ", outDir)
