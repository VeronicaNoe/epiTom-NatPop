suppressPackageStartupMessages(library(data.table))
options(scipen = 0)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) {
  stop("Usage: Rscript fig5i_OLP_state_frequencies.R [meth_table] [sample_table] [out_dir]")
}
meth_file <- if (length(args) >= 1L) args[1L] else
  "/mnt/disk2/vibanez/10_data-analysis/Fig5/results/02.0_merged-mLevels-mStatus-metaLevels.tsv"
sample_file <- if (length(args) >= 2L) args[2L] else
  "/mnt/disk2/vibanez/10_data-analysis/Fig1/00_sampleTab.tsv"
out_dir <- if (length(args) >= 3L) args[3L] else
  "/mnt/disk2/vibanez/10_data-analysis/Fig5/results/Tm122_OLP_state_frequency"
metabolite_id <- "Tm122"
target_dmr <- "ch12_C-DMR_62312001"
run_pairwise_tests <- TRUE
if (!file.exists(meth_file)) stop("Missing file: ", meth_file)
if (!file.exists(sample_file)) stop("Missing file: ", sample_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Cannot create output directory: ", out_dir)
write_tsv <- function(x, name) {
  fwrite(x, file.path(out_dir, name), sep = "\t", quote = FALSE, na = "NA")
}
clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  x
}
raw <- fread(meth_file, sep = "\t", header = TRUE, fill = TRUE,
             na.strings = c("NA", "NaN", ""),
             select = c("ID", "Accessions", "methStatus", "metabolite"))
x <- raw[ID == target_dmr & metabolite == metabolite_id,
         .(Accessions, methStatus)]
rm(raw)
if (!nrow(x)) stop("No rows found for ", metabolite_id, " and ", target_dmr)
x[, Accessions := clean_text(Accessions)]
x[, methStatus := clean_text(methStatus)]
if (anyNA(x$Accessions)) stop("Some target rows lack accession identifiers.")
cat("Target rows before removing repeated annotations:", nrow(x), "\n")
print(x[, .N, by = methStatus])
unknown <- setdiff(unique(x$methStatus[!is.na(x$methStatus)]), c("0", "2"))
if (length(unknown)) {
  write_tsv(x, "CHECK_target_rows_unknown_states.tsv")
  stop("Unexpected methStatus codes: ", paste(unknown, collapse = ", "),
       ". Their meaning must be established before choosing a denominator.")
}
conflicts <- x[!is.na(methStatus),
               .(n_states = uniqueN(methStatus)), by = Accessions][n_states > 1L]
if (nrow(conflicts)) {
  write_tsv(x[Accessions %chin% conflicts$Accessions], "CHECK_conflicting_states.tsv")
  stop("Some accessions have conflicting state calls. See CHECK_conflicting_states.tsv.")
}
states <- x[, .(methStatus = {
  v <- unique(methStatus[!is.na(methStatus)])
  if (length(v)) v[1L] else NA_character_
}), by = Accessions]
samples <- fread(sample_file, sep = "\t", header = TRUE, quote = "\"",
                 fill = TRUE, na.strings = c("NA", "NaN", ""),
                 select = c("AccName", "Group"))
samples[, AccName := clean_text(AccName)]
samples[, Group := clean_text(Group)]
samples <- unique(samples[!is.na(AccName) & !is.na(Group),
                          .(Accessions = AccName, Group)])
group_order <- c("WILD", "PIM", "SLC", "SLL")
groups <- c(intersect(group_order, unique(samples$Group)),
            sort(setdiff(unique(samples$Group), group_order)))
sample_map <- samples[Accessions %chin% states$Accessions]
conflicting_groups <- sample_map[, .(n_groups = uniqueN(Group)),
                                 by = Accessions][n_groups > 1L]
if (nrow(conflicting_groups)) {
  write_tsv(sample_map[Accessions %chin% conflicting_groups$Accessions],
            "CHECK_conflicting_groups.tsv")
  stop("Conflicting population labels. See CHECK_conflicting_groups.tsv.")
}
unmatched <- states[!Accessions %chin% sample_map$Accessions]
if (nrow(unmatched)) {
  write_tsv(unmatched, "CHECK_unmatched_accessions.tsv")
}
acc <- merge(states, sample_map, by = "Accessions", all.x = TRUE)
acc[, state := fifelse(is.na(methStatus), NA_character_,
                       fifelse(methStatus == "0", "UM", "M"))]
acc[, `:=`(ID = target_dmr, metabolite = metabolite_id)]
setcolorder(acc, c("ID", "metabolite", "Accessions", "Group", "methStatus", "state"))
setorder(acc, Group, Accessions)
write_tsv(acc, "01_accession_states.tsv")
summary <- acc[, .(
  n_present = .N,
  n_missing = sum(is.na(state)),
  N = sum(!is.na(state)),
  n_UM = sum(state == "UM", na.rm = TRUE),
  n_M = sum(state == "M", na.rm = TRUE)
), by = Group]
summary <- merge(data.table(Group = groups), summary, by = "Group", all.x = TRUE)
count_cols <- c("n_present", "n_missing", "N", "n_UM", "n_M")
for (nm in count_cols) set(summary, which(is.na(summary[[nm]])), nm, 0L)
summary <- summary[match(groups, Group)]
summary[, `:=`(
  frequency_UM = fifelse(N > 0, n_UM / N, NA_real_),
  frequency_M = fifelse(N > 0, n_M / N, NA_real_)
)]
stopifnot(all(summary$N == summary$n_UM + summary$n_M))
stopifnot(all(summary$n_present == summary$N + summary$n_missing))
write_tsv(summary, "02_population_state_frequencies.tsv")
cat("\nFrequencies by population (fractions, not percentages):\n")
print(summary)
if (run_pairwise_tests) {
  available <- summary[N > 0, Group]
  if (length(available) >= 2L) {
    pairs <- combn(available, 2L, simplify = FALSE)
    tests <- rbindlist(lapply(pairs, function(pair) {
      a <- summary[Group == pair[1L]]
      b <- summary[Group == pair[2L]]
      ct <- matrix(c(a$n_M, a$n_UM, b$n_M, b$n_UM), nrow = 2L, byrow = TRUE,
                   dimnames = list(pair, c("M", "UM")))
      monomorphic <- any(colSums(ct) == 0L)
      if (monomorphic) {
        p <- 1
        odds <- lower <- upper <- NA_real_
      } else {
        fit <- fisher.test(ct, alternative = "two.sided", conf.level = 0.95)
        p <- fit$p.value
        odds <- unname(fit$estimate)
        lower <- fit$conf.int[1L]
        upper <- fit$conf.int[2L]
      }
      data.table(group1 = pair[1L], group2 = pair[2L],
                 M1 = a$n_M, UM1 = a$n_UM, N1 = a$N,
                 M2 = b$n_M, UM2 = b$n_UM, N2 = b$N,
                 odds_ratio_M_group1_vs_group2 = odds,
                 lower95_OR = lower, upper95_OR = upper, P = p,
                 note = if (monomorphic) "One state absent in both groups; OR undefined" else "")
    }))
    tests[, P_BH := p.adjust(P, method = "BH")]
    write_tsv(tests, "03_population_pairwise_Fisher_BH.tsv")
    cat("\nTwo-sided Fisher tests; BH adjustment across all", nrow(tests), "population pairs:\n")
    print(tests)
  } else {
    message("Fewer than two populations have valid state calls: no pairwise tests performed.")
  }
}
message("Saved results in: ", out_dir,
        "\nThese frequencies describe the available Tm122/OLP accession subset, not necessarily the full panel.")
