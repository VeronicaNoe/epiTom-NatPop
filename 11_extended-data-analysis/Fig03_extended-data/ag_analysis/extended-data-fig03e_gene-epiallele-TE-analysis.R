#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})


## ==========================================================================
## 0. Paths and settings
## ==========================================================================

inputDir <- "/mnt/disk2/vibanez/10_data-analysis/Fig4/results"
outDir <- "/mnt/disk2/vibanez/11_extended-data-analysis/Fig03_extended-data"

allDataFile <- file.path(
  inputDir,
  "02.07_epiallele_sample.tsv"
)

closestTEFile <- file.path(
  outDir,
  "extended-data-fig03e_all-genes_closest-TE.tsv"
)

classes <- c(
  "UM",
  "gbM",
  "teM"
)

allClassLevels <- c(
  "Other",
  classes
)

binWidth <- 500L
maxDistance <- 3000L

classColours <- c(
  "UM"  = "grey50",
  "gbM" = "#ffca7b",
  "teM" = "#820a86"
)

barOffset <- 0.21
barWidth <- 0.38

dir.create(
  outDir,
  recursive = TRUE,
  showWarnings = FALSE
)

## ==========================================================================
## 1. Helper functions
## ==========================================================================

format_p <- function(p) {
  if (is.na(p)) {
    return("NA")
  }

  if (p < 1e-4) {
    return("<1e-4")
  }

  if (p < 0.001) {
    return(
      format(
        p,
        scientific = TRUE,
        digits = 2
      )
    )
  }

  sprintf(
    "%.4f",
    p
  )
}


significance_label <- function(p) {
  fifelse(
    is.na(p),
    "",
    fifelse(
      p < 0.001,
      "***",
      fifelse(
        p < 0.01,
        "**",
        fifelse(
          p < 0.05,
          "*",
          ""
        )
      )
    )
  )
}


require_columns <- function(
    x,
    required,
    objectName
) {
  missingColumns <- setdiff(
    required,
    names(x)
  )

  if (length(missingColumns) > 0L) {
    stop(
      objectName,
      " is missing required columns: ",
      paste(
        missingColumns,
        collapse = ", "
      )
    )
  }
}


make_distance_definition <- function(
    maxDistance,
    binWidth
) {
  upstreamKeys <- as.character(
    seq(
      -maxDistance,
      -binWidth,
      by = binWidth
    )
  )

  downstreamKeys <- as.character(
    seq(
      binWidth,
      maxDistance,
      by = binWidth
    )
  )

  keys <- c(
    upstreamKeys,
    "Overlap",
    downstreamKeys,
    "Outside"
  )

  labels <- c(
    paste0(
      seq(
        -maxDistance,
        -binWidth,
        by = binWidth
      ) / 1000
    ),
    "Overlap",
    paste0(
      seq(
        binWidth,
        maxDistance,
        by = binWidth
      ) / 1000
    ),
    ">3 kb / no TE"
  )

  names(
    labels
  ) <- keys

  list(
    keys = keys,
    labels = labels
  )
}


safe_fisher <- function(
    targetInBin,
    targetOutsideBin,
    referenceInBin,
    referenceOutsideBin
) {
  testTable <- matrix(
    c(
      targetInBin,
      targetOutsideBin,
      referenceInBin,
      referenceOutsideBin
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c(
        "Target",
        "Reference"
      ),
      c(
        "In_bin",
        "Outside_bin"
      )
    )
  )

  test <- fisher.test(
    testTable,
    alternative = "two.sided"
  )

  data.table(
    oddsRatio = unname(
      test$estimate
    ),
    lower95OR = unname(
      test$conf.int[1]
    ),
    upper95OR = unname(
      test$conf.int[2]
    ),
    pValue = test$p.value
  )
}


build_occurrence_test <- function(
    geneData,
    classes,
    activeDistanceKeys,
    activeDistanceLabels,
    comparison = c(
      "Other",
      "Genome"
    )
) {
  comparison <- match.arg(
    comparison
  )

  resultList <- vector(
    "list",
    length(
      classes
    ) *
      length(
        activeDistanceKeys
      )
  )

  resultIndex <- 1L

  for (
    className in classes
  ) {
    if (
      comparison == "Other"
    ) {
      analysisData <- copy(
        geneData[
          as.character(
            geneClass
          ) %chin% c(
            className,
            "Other"
          )
        ]
      )

      analysisData[
        ,
        target := as.character(
          geneClass
        ) == className
      ]

      referenceLabel <- "Other"
    } else {
      analysisData <- copy(
        geneData
      )

      analysisData[
        ,
        target := as.character(
          geneClass
        ) == className
      ]

      referenceLabel <- "All non-class genes"
    }

    targetTotal <- analysisData[
      target == TRUE,
      .N
    ]

    referenceTotal <- analysisData[
      target == FALSE,
      .N
    ]

    universeTotal <- targetTotal +
      referenceTotal

    for (
      distanceName in activeDistanceKeys
    ) {
      targetInBin <- analysisData[
        target == TRUE &
          distanceKey == distanceName,
        .N
      ]

      referenceInBin <- analysisData[
        target == FALSE &
          distanceKey == distanceName,
        .N
      ]

      targetOutsideBin <- targetTotal -
        targetInBin

      referenceOutsideBin <- referenceTotal -
        referenceInBin

      binTotal <- targetInBin +
        referenceInBin

      expectedMean <- targetTotal *
        binTotal /
        universeTotal

      # Exact distribution of the number of target genes obtained when a
      # random set of targetTotal genes is drawn from this universe.
      expectedLower95 <- qhyper(
        0.025,
        m = binTotal,
        n = universeTotal -
          binTotal,
        k = targetTotal
      )

      expectedUpper95 <- qhyper(
        0.975,
        m = binTotal,
        n = universeTotal -
          binTotal,
        k = targetTotal
      )

      fisherResult <- safe_fisher(
        targetInBin,
        targetOutsideBin,
        referenceInBin,
        referenceOutsideBin
      )

      resultList[[resultIndex]] <- data.table(
        comparison = comparison,
        reference = referenceLabel,
        geneClass = className,
        distanceKey = distanceName,
        targetTotal = targetTotal,
        referenceTotal = referenceTotal,
        universeTotal = universeTotal,
        targetInBin = targetInBin,
        referenceInBin = referenceInBin,
        observedCount = targetInBin,
        expectedMean = expectedMean,
        expectedLower95 = expectedLower95,
        expectedUpper95 = expectedUpper95,
        oddsRatio = fisherResult$oddsRatio,
        lower95OR = fisherResult$lower95OR,
        upper95OR = fisherResult$upper95OR,
        pValue = fisherResult$pValue
      )

      resultIndex <- resultIndex +
        1L
    }
  }

  result <- rbindlist(
    resultList
  )

  result[
    ,
    adjustedP := p.adjust(
      pValue,
      method = "BH"
    ),
    by = geneClass
  ]

  result[
    ,
    `:=`(
      geneClass = factor(
        geneClass,
        levels = classes
      ),
      distanceCategory = factor(
        distanceKey,
        levels = activeDistanceKeys,
        labels = activeDistanceLabels
      ),
      xIndex = match(
        distanceKey,
        activeDistanceKeys
      ),
      significance = significance_label(
        adjustedP
      ),
      direction = fcase(
        oddsRatio > 1,
        "Enriched",

        oddsRatio < 1,
        "Depleted",

        default = "No difference"
      )
    )
  ]

  result[
    ,
    panelMaximum := max(
      c(
        observedCount,
        expectedUpper95
      ),
      na.rm = TRUE
    ),
    by = geneClass
  ]

  result[
    ,
    `:=`(
      observedX = xIndex -
        barOffset,
      expectedX = xIndex +
        barOffset,
      observedLabel = comma(
        observedCount
      ),
      expectedLabel = comma(
        round(
          expectedMean
        )
      ),
      starY = pmax(
        observedCount,
        expectedUpper95
      ) +
        0.045 *
        panelMaximum
    )
  ]

  setorder(
    result,
    geneClass,
    xIndex
  )

  result
}


plot_occurrence_test <- function(
    result,
    captionText,
    outputFile,
    overlapIndex,
    activeDistanceKeys,
    activeDistanceLabels
) {
  plotObject <- ggplot(
    result
  ) +
    geom_col(
      aes(
        x = observedX,
        y = observedCount,
        fill = geneClass,
        colour = geneClass
      ),
      width = barWidth,
      linewidth = 0.35,
      show.legend = FALSE
    ) +
    geom_col(
      aes(
        x = expectedX,
        y = expectedMean,
        colour = geneClass
      ),
      fill = "white",
      width = barWidth,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_errorbar(
      aes(
        x = expectedX,
        ymin = expectedLower95,
        ymax = expectedUpper95,
        colour = geneClass
      ),
      width = 0.10,
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    geom_text(
      aes(
        x = observedX,
        y = observedCount,
        label = observedLabel
      ),
      vjust = -0.35,
      size = 2.25
    ) +
    geom_text(
      aes(
        x = expectedX,
        y = expectedUpper95,
        label = expectedLabel
      ),
      vjust = -0.35,
      size = 2.25
    ) +
    geom_text(
      aes(
        x = observedX,
        y = starY,
        label = significance
      ),
      size = 3.4,
      fontface = "bold"
    ) +
    geom_vline(
      xintercept = overlapIndex,
      linetype = "dotted",
      linewidth = 0.4
    ) +
    facet_wrap(
      ~geneClass,
      ncol = 1,
      scales = "free_y"
    ) +
    scale_fill_manual(
      values = classColours,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = classColours,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(
        activeDistanceKeys
      ),
      labels = activeDistanceLabels,
      expand = expansion(
        mult = c(
          0.02,
          0.03
        )
      )
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(
        mult = c(
          0,
          0.15
        )
      )
    ) +
    labs(
      x = "Closest-TE position relative to the gene",
      y = "Number of genes",
      caption = captionText
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.caption = element_text(
        size = 8,
        hjust = 0
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold"
      )
    )

  print(
    plotObject
  )

  ggsave(
    outputFile,
    plotObject,
    width = 25,
    height = 25,
    units = "cm"
  )

  invisible(
    plotObject
  )
}


## ==========================================================================
## 2. Read and validate data
## ==========================================================================

if (!file.exists(allDataFile)) {
  stop(
    "Missing input file: ",
    allDataFile
  )
}

if (!file.exists(closestTEFile)) {
  stop(
    "Missing closest-TE file: ",
    closestTEFile
  )
}

allData <- fread(
  allDataFile,
  sep = "\t",
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = c(
    "NA",
    ""
  ),
  nThread = 10
)

allGeneTE <- fread(
  closestTEFile,
  sep = "\t",
  data.table = TRUE,
  fill = TRUE,
  check.names = FALSE,
  na.strings = c(
    "NA",
    ""
  ),
  nThread = 10
)

setDT(
  allData
)

setDT(
  allGeneTE
)

require_columns(
  allData,
  c(
    "geneName",
    "geneEpiallele"
  ),
  "allData"
)

require_columns(
  allGeneTE,
  c(
    "geneName",
    "geneStrand",
    "teChr",
    "teStart",
    "teEnd",
    "signedDistance"
  ),
  "allGeneTE"
)

allGeneTE[
  ,
  `:=`(
    teStart = suppressWarnings(
      as.numeric(
        teStart
      )
    ),
    teEnd = suppressWarnings(
      as.numeric(
        teEnd
      )
    ),
    signedDistance = suppressWarnings(
      as.numeric(
        signedDistance
      )
    )
  )
]
## ==========================================================================
## 3. Assign all genes to Other, UM, gbM or teM
## ==========================================================================

classMap <- unique(
  allData[
    geneEpiallele %chin% classes &
      !is.na(
        geneName
      ),
    .(
      geneName,
      geneEpiallele
    )
  ]
)

classConflicts <- classMap[
  ,
  .(
    nClasses = uniqueN(
      geneEpiallele
    )
  ),
  by = geneName
][
  nClasses > 1L
]

if (
  nrow(
    classConflicts
  ) > 0L
) {
  print(
    classConflicts
  )


}

duplicatedClosestGenes <- allGeneTE[
  ,
  .N,
  by = geneName
][
  N > 1L
]

if (
  nrow(
    duplicatedClosestGenes
  ) > 0L
) {
  print(
    duplicatedClosestGenes
  )

  stop(
    "The closest-TE table contains more than one row per gene."
  )
}

geneData <- merge(
  allGeneTE,
  classMap,
  by = "geneName",
  all.x = TRUE
)

geneData[
  is.na(
    geneEpiallele
  ),
  geneEpiallele := "Other"
]

geneData[
  ,
  geneClass := factor(
    geneEpiallele,
    levels = allClassLevels
  )
]

geneData[
  ,
  geneEpiallele := NULL
]


## ==========================================================================
## 4. Define closest-TE categories
## ==========================================================================

distanceDefinition <- make_distance_definition(
  maxDistance = maxDistance,
  binWidth = binWidth
)

allDistanceKeys <- distanceDefinition$keys
allDistanceLabels <- distanceDefinition$labels

geneData[
  ,
  noTEHit := is.na(
    teChr
  ) |
    teChr == "." |
    is.na(
      teStart
    ) |
    teStart < 0
]

geneData[
  ,
  distanceKey := fcase(
    noTEHit,
    "Outside",

    abs(
      signedDistance
    ) > maxDistance,
    "Outside",

    signedDistance == 0,
    "Overlap",

    default = as.character(
      sign(
        signedDistance
      ) *
        binWidth *
        ceiling(
          abs(
            signedDistance
          ) / binWidth
        )
    )
  )
]

invalidDistance <- geneData[
  !distanceKey %chin% allDistanceKeys
]

if (
  nrow(
    invalidDistance
  ) > 0L
) {
  print(
    invalidDistance
  )


}

activeDistanceKeys <- allDistanceKeys[
  allDistanceKeys %chin% unique(
    geneData$distanceKey
  )
]

activeDistanceLabels <- unname(
  allDistanceLabels[
    activeDistanceKeys
  ]
)

geneData[
  ,
  `:=`(
    distanceCategory = factor(
      distanceKey,
      levels = activeDistanceKeys,
      labels = activeDistanceLabels
    ),
    xIndex = match(
      distanceKey,
      activeDistanceKeys
    )
  )
]

overlapIndex <- match(
  "Overlap",
  activeDistanceKeys
)


## ==========================================================================
## 5. Diagnostics
## ==========================================================================

cat(
  "\n============================================================\n"
)

cat(
  "GENE UNIVERSE\n"
)

cat(
  "============================================================\n"
)

print(
  geneData[
    ,
    .(
      nGenes = .N,
      percentage = 100 *
        .N /
        nrow(
          geneData
        )
    ),
    by = geneClass
  ]
)

cat(
  "\nClosest-TE categories:\n"
)

print(
  geneData[
    ,
    .(
      nGenes = .N,
      percentage = 100 *
        .N /
        nrow(
          geneData
        )
    ),
    by = distanceCategory
  ][
    order(
      distanceCategory
    )
  ]
)


###############################################################################
## Genome-wide enrichment of each epiallele class
###############################################################################
##
## Each class is compared with all genes that do not belong to that class.
##
## The expected bar is the exact expectation for a random gene set drawn from
## the complete annotated-gene universe with the same size as the class.
##
## This answers:
##
##   Are there more or fewer class genes in this TE-distance bin than expected
##   from the genome-wide gene distribution?
##
## Only per-bin tests are performed; there is no global/omnibus test.

genomeEnrichment <- build_occurrence_test(
  geneData = geneData,
  classes = classes,
  activeDistanceKeys = activeDistanceKeys,
  activeDistanceLabels = activeDistanceLabels,
  comparison = "Genome"
)

fwrite(
  genomeEnrichment,
  file.path(
    outDir,
    "extended-data-fig03e_genome-enrichment_per-bin-tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "\n============================================================\n"
)

cat(
  "GENOME-WIDE EPIALLELE-CLASS ENRICHMENT\n"
)

cat(
  "============================================================\n"
)

print(
  genomeEnrichment[
    ,
    .(
      geneClass,
      distanceCategory,
      observedCount,
      expectedMean,
      expectedLower95,
      expectedUpper95,
      oddsRatio,
      lower95OR,
      upper95OR,
      direction,
      pValue,
      adjustedP,
      significance
    )
  ]
)

pGenomeEnrichment <- plot_occurrence_test(
  result = genomeEnrichment,
  captionText = paste0(
    "Filled bars: observed counts. Open bars: expected counts for a random ",
    "gene set of the same size drawn from the complete annotated-gene universe. ",
    "Error bars are exact hypergeometric 95% intervals. ",
    "Per-bin Fisher exact tests comparing the class with all non-class genes; ",
    "BH correction across bins within each class. ",
    "* adjusted P < 0.05; ** < 0.01; *** < 0.001."
  ),
  outputFile = file.path(
    outDir,
    "extended-data-fig03e_genome-enrichment_observed-expected-counts.pdf"
  ),
  overlapIndex = overlapIndex,
  activeDistanceKeys = activeDistanceKeys,
  activeDistanceLabels = activeDistanceLabels
)
