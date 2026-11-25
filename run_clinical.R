#!/usr/bin/env Rscript

# Some Windows command shells start R in the C locale even for Unicode paths.
if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
}
args <- commandArgs(trailingOnly = TRUE)
usage <- paste("Rscript run_clinical.R <normalized_nCounter.xlsx> <patient_annotation.xlsx>",
               "<new_output_directory> <sample_metadata.csv|Yamauchi_Supplementary_Data_3.xlsx>")
if (length(args) == 1L && args %in% c("--help", "-h")) {
  cat(usage, "\n")
  quit(status = 0)
}
if (length(args) != 4L) stop(usage, call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
inputs <- normalizePath(args[c(1, 2, 4)], mustWork = TRUE)
if (any(dir.exists(inputs))) stop("Each input must be a file.")
packages <- c("readxl", "ggplot2", "pheatmap", "survival")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
clinical_sheets <- readxl::excel_sheets(inputs[2])
clinical_sheet <- if ("Patient annotation" %in% clinical_sheets) "Patient annotation" else "Fig1b_annotation"
if (!clinical_sheet %in% clinical_sheets) stop("Missing Patient annotation or Fig1b_annotation sheet.")
metadata_format <- tolower(tools::file_ext(inputs[3]))
if (!metadata_format %in% c("csv", "xlsx")) stop("Expected a metadata CSV or the source-data XLSX.")
out <- args[3]
if (file.exists(out)) stop("Choose a new output directory.")
if (!dir.exists(dirname(out)) && !dir.create(dirname(out), recursive = TRUE)) stop("Cannot create output parent.")
stage <- tempfile(pattern = ".clinical-run-", tmpdir = dirname(out))
if (!dir.create(stage)) stop("Cannot prepare output directory.")
message("Writing to ", stage, " until all analyses finish.")
metadata <- inputs[3]
if (metadata_format == "xlsx") {
  # Yamauchi et al. identify the eight week-4 pairs in the source data for their Figure 3b.
  rows <- suppressMessages(readxl::read_xlsx(inputs[1], sheet = "NC_9785_N1_Normalization",
    col_names = FALSE, col_types = "text", n_max = 3))
  rcc <- grepl("\\.RCC$", unlist(rows[1, ], use.names = FALSE), ignore.case = TRUE)
  ids <- trimws(unlist(rows[3, rcc], use.names = FALSE))
  header <- suppressMessages(readxl::read_xlsx(metadata, sheet = "Fig3b", skip = 1, n_max = 1,
    col_names = FALSE, col_types = "text"))
  labels <- trimws(unlist(header, use.names = FALSE))
  week4 <- labels[grepl("^FL[0-9]+_POST$", labels)]
  week4_baseline <- labels[grepl("^FL[0-9]+$", labels)]
  if (!all(c("mean_Pre", "mean_4W") %in% labels) || length(week4) != 8L ||
      anyDuplicated(week4) || !setequal(week4_baseline, sub("_POST$", "", week4)) ||
      !all(week4 %in% ids)) stop("Unexpected week-4 sample identifiers in Fig3b.")
  if (length(ids) != 59L || anyDuplicated(ids) || anyNA(ids) ||
      !all(grepl("^FL[0-9]+(_POST)?$", ids))) stop("Expected 59 unique Yamauchi cohort sample IDs.")
  followup <- endsWith(ids, "_POST")
  if (sum(followup) != 20L || !all(sub("_POST$", "", ids[followup]) %in% ids[!followup])) {
    stop("Expected 20 follow-up biopsies with matching baseline IDs.")
  }
  annotation <- data.frame(sample_id = ids, patient_id = sub("_POST$", "", ids),
    timepoint = ifelse(!followup, "baseline", ifelse(ids %in% week4, "week4", "PD")))
  metadata <- file.path(stage, "sample_metadata.csv")
  utils::write.csv(annotation, metadata, row.names = FALSE)
}
run <- function(script, ...) {
  implementation <- new.env(parent = globalenv())
  sys.source(file.path(here, "R", "Human", script), envir = implementation)
  implementation$main(c(...))
}
scores <- file.path(stage, "step1", "signature_scores_and_z_by_sample.csv")
deltas <- file.path(stage, "step2", "PD_delta_long_by_patient_signature.csv")
statistics <- file.path(stage, "step2b", "S1B_STable_PD_delta_stats_ALLsignatures.csv")
run("Step1_Human_nCounter_signature_scores.R", inputs[1], metadata, file.path(stage, "step1"))
run("Step2_Human_PD_delta_heatmaps.R", inputs[2], scores, file.path(stage, "step2"), clinical_sheet)
run("Step2b_Human_PD_signature_statistics.R", deltas, file.path(stage, "step2b"))
run("Step3_Human_PD_signature_ranking.R", statistics, file.path(stage, "step3"))
run("Paired_Human_nCounter_panels.R", inputs[1], metadata, inputs[2], file.path(stage, "paired"), clinical_sheet)
run("Human_Fig1C_PFS.R", scores, inputs[2], file.path(stage, "pfs"), clinical_sheet)
if (file.exists(out) || !file.rename(stage, out)) stop("Completed outputs remain in ", stage)
message("Clinical outputs: ", normalizePath(out, winslash = "/"))
