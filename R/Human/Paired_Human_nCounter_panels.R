#!/usr/bin/env Rscript

# Fig. 1B and S1C-F: paired-sample, gene-level standardization.
# The 20 patients with paired biopsies define the 40-sample reference. Only
# the 12 patients with progression biopsies enter the plotted comparisons.

read_expression <- function(path, genes) {
  raw <- suppressMessages(readxl::read_xlsx(path, sheet = "NC_9785_N1_Normalization",
                                            col_names = FALSE, col_types = "text"))
  header <- trimws(as.character(unlist(raw[1, ], use.names = FALSE)))
  sample_cols <- which(!is.na(header) & grepl("\\.RCC$", header, ignore.case = TRUE))
  if (length(sample_cols) != 59L || header[1] != "Probe Name") stop("Unexpected nCounter header.")
  sample_ids <- trimws(as.character(unlist(raw[3, sample_cols], use.names = FALSE)))
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Sample IDs in row 3 must be unique and nonempty.")
  }
  probe <- trimws(as.character(raw[[1]]))
  if (!all(genes %in% probe) || anyDuplicated(probe[probe %in% genes])) {
    stop("Each marker must have one expression row.")
  }
  x <- as.matrix(raw[match(genes, probe), sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(x) <- "double")
  if (any(!is.finite(x)) || any(x < 0)) stop("Expected nonnegative, finite marker expression.")
  dimnames(x) <- list(genes, sample_ids)
  x
}

read_pairs <- function(path, samples) {
  meta <- utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
  keys <- c("patient_id", "sample_id", "timepoint")
  if (!all(keys %in% names(meta))) stop("Missing sample metadata columns.")
  meta <- meta[keys]
  meta[] <- lapply(meta, trimws)
  if (anyNA(meta) || any(as.matrix(meta) == "") || anyDuplicated(meta$sample_id) ||
      anyDuplicated(meta[c("patient_id", "timepoint")]) || !setequal(samples, meta$sample_id) ||
      !all(meta$timepoint %in% c("baseline", "week4", "PD"))) stop("Invalid sample metadata.")
  follow <- meta[meta$timepoint != "baseline", ]
  base <- meta[meta$timepoint == "baseline", ]
  if (nrow(base) != 39L || anyDuplicated(follow$patient_id) || sum(follow$timepoint == "week4") != 8L ||
      sum(follow$timepoint == "PD") != 12L || !all(follow$patient_id %in% base$patient_id)) {
    stop("Expected 20 distinct patients with baseline and either week-4 or progression biopsies.")
  }
  follow <- follow[order(follow$patient_id), ]
  data.frame(patient_id = follow$patient_id, baseline_sample = base$sample_id[match(follow$patient_id, base$patient_id)],
             followup_sample = follow$sample_id, timepoint = follow$timepoint)
}

read_response_groups <- function(path, ids, sheet) {
  ann <- as.data.frame(readxl::read_xlsx(path, sheet = sheet,
    skip = as.integer(sheet == "Fig1b_annotation")))
  if (!all(c("FL_ID", "mRECIST") %in% names(ann))) stop("Missing FL_ID or mRECIST.")
  patient <- trimws(as.character(ann$FL_ID))
  if (anyDuplicated(patient[patient %in% ids])) stop("Duplicate paired patient in clinical annotation.")
  response <- toupper(trimws(as.character(ann$mRECIST[match(ids, patient)])))
  if (!all(response %in% c("PD", "CR", "PR", "SD"))) stop("Missing or unknown mRECIST response.")
  group <- ifelse(response == "PD", "PD", "non-PD")
  if (sum(group == "PD") != 3L || sum(group == "non-PD") != 9L) stop("Unexpected response group sizes.")
  group
}

plot_pairs <- function(coords, tests, labels, path) {
  grDevices::pdf(path, width = 7, height = 4.7, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  for (panel in names(labels)) {
    dat <- coords[coords$panel == panel, ]
    graphics::par(mfrow = c(1, 2), mar = c(3.2, 4.1, 2.7, 1), oma = c(0, 0, 1.5, 0))
    limits <- range(c(dat$baseline, dat$progression))
    for (g in c("non-PD", "PD")) {
      d <- dat[dat$group == g, ]
      p <- tests$p_value[tests$panel == panel & tests$group == g]
      graphics::plot(NA_real_, xlim = c(0.8, 2.2), ylim = limits, xaxt = "n", xlab = "",
                     ylab = labels[[panel]], main = sprintf("%s (n = %d), P = %.4f", g, nrow(d), p))
      graphics::axis(1, at = c(1, 2), labels = c("Baseline", "Progression"))
      graphics::segments(1, d$baseline, 2, d$progression, col = "grey55")
      graphics::points(rep(1, nrow(d)), d$baseline, pch = 16)
      graphics::points(rep(2, nrow(d)), d$progression, pch = 16)
    }
    graphics::mtext(paste("Fig.", panel), side = 3, outer = TRUE, font = 2)
  }
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Usage: Rscript Paired_Human_nCounter_panels.R <nCounter.xlsx>",
                 "<sample_metadata.csv> <patient_annotation.xlsx> <new_output_dir> [clinical_sheet]")
  if (length(args) == 1L && args[1] %in% c("-h", "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
  if (!length(args) %in% 4:5) stop(usage, call. = FALSE)
  sheet <- if (length(args) == 5L) args[5] else "Patient annotation"
  out <- args[4]
  if (file.exists(out)) stop("Choose a new output directory.")
  sets <- list(Neutrophil = c("CXCR2", "CSF3R"), Macrophage = c("CD68", "CSF1R"), CD8_T = c("CD8A", "CD8B"))
  x <- read_expression(args[1], unique(unlist(sets, use.names = FALSE)))
  pairs <- read_pairs(args[2], colnames(x))
  x <- log2(1 + x[, c(pairs$baseline_sample, pairs$followup_sample), drop = FALSE])
  mu <- rowMeans(x)
  centered <- sweep(x, 1, mu, "-")
  pop_sd <- sqrt(rowMeans(centered^2))
  if (any(pop_sd <= 0)) stop("Cannot standardize a constant marker.")
  z <- sweep(centered, 1, pop_sd, "/")
  sig <- vapply(sets, function(genes) colMeans(z[genes, , drop = FALSE]), numeric(ncol(z)))
  scores <- cbind("1B" = sig[, "Neutrophil"] - sig[, "CD8_T"], "S1C" = sig[, "Neutrophil"],
                  "S1D" = sig[, "Macrophage"], "S1E" = sig[, "CD8_T"], "S1F" = sig[, "Macrophage"] - sig[, "CD8_T"])
  rownames(scores) <- colnames(z)
  shown <- pairs[pairs$timepoint == "PD", ]
  group <- read_response_groups(args[3], shown$patient_id, sheet)
  coords <- do.call(rbind, lapply(colnames(scores), function(panel) {
    baseline <- unname(scores[shown$baseline_sample, panel])
    progression <- unname(scores[shown$followup_sample, panel])
    data.frame(panel, patient_id = shown$patient_id, group, baseline, progression, delta = progression - baseline)
  }))
  tests <- do.call(rbind, lapply(colnames(scores), function(panel) {
    do.call(rbind, lapply(c("non-PD", "PD"), function(g) {
      d <- coords[coords$panel == panel & coords$group == g, ]
      t <- stats::t.test(d$progression, d$baseline, paired = TRUE, alternative = "two.sided")
      data.frame(panel, group = g, n = nrow(d), mean_delta = mean(d$delta),
                 t = unname(t$statistic), df = unname(t$parameter), p_value = t$p.value)
    }))
  }))
  parent <- dirname(out)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) stop("Cannot create output parent.")
  stage <- tempfile(".paired_panels_", tmpdir = parent)
  if (!dir.create(stage)) stop("Cannot create staging directory.")
  reference <- data.frame(marker = names(mu), reference_n = ncol(z), mean_log2 = unname(mu), population_sd = unname(pop_sd))
  utils::write.csv(reference, file.path(stage, "paired_marker_reference.csv"), row.names = FALSE)
  utils::write.csv(coords, file.path(stage, "paired_panel_coordinates.csv"), row.names = FALSE)
  utils::write.csv(tests, file.path(stage, "paired_panel_tests.csv"), row.names = FALSE)
  labels <- c("1B" = "Neutrophil Z - CD8 T Z", "S1C" = "Neutrophil Z", "S1D" = "Macrophage Z",
              "S1E" = "CD8 T Z", "S1F" = "Macrophage Z - CD8 T Z")
  plot_pairs(coords, tests, labels, file.path(stage, "paired_panel_plots.pdf"))
  writeLines(capture.output(sessionInfo()), file.path(stage, "sessionInfo.txt"))
  if (file.exists(out) || !file.rename(stage, out)) stop("Could not finalize output; staged files are in ", stage)
  message("Wrote five panels and ten paired tests to ", out)
}

if (sys.nframe() == 0L) main()
