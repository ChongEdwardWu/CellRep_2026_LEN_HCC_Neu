#!/usr/bin/env Rscript

# Reconstruct the published event curve using the direction of the neutrophil
# signature change between baseline and the observed on-treatment biopsy.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Usage: Rscript Human_Fig1C_PFS.R <signature_scores_and_z_by_sample.csv>",
                 "<patient_annotation.xlsx> <new_output_directory> [clinical_sheet]")
  if (length(args) == 1L && args %in% c("-h", "--help")) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  if (!length(args) %in% 3:4) stop(usage, call. = FALSE)
  sheet <- if (length(args) == 4L) args[4] else "Patient annotation"
  scores <- utils::read.csv(args[1], stringsAsFactors = FALSE)
  required <- c("patient_id", "timepoint", "Neutrophil_z")
  if (!all(required %in% names(scores)) || anyNA(scores[required]) ||
      anyDuplicated(scores[c("patient_id", "timepoint")]) || !all(is.finite(scores$Neutrophil_z))) {
    stop("Invalid Step 1 neutrophil signature table.")
  }
  baseline <- scores[scores$timepoint == "baseline", required]
  post <- scores[scores$timepoint %in% c("week4", "PD"), required]
  if (!nrow(post) || anyDuplicated(post$patient_id) ||
      !all(post$patient_id %in% baseline$patient_id)) stop("Expected one paired follow-up per patient.")
  pre <- baseline[match(post$patient_id, baseline$patient_id), ]
  delta <- post$Neutrophil_z - pre$Neutrophil_z
  if (any(delta == 0)) stop("A zero change falls outside the two directions defined in Figure 1C.")
  clinical <- as.data.frame(readxl::read_xlsx(args[2], sheet = sheet,
    skip = as.integer(sheet == "Fig1b_annotation")))
  fields <- c("FL_ID", "PFS (month)", "censor")
  if (!all(fields %in% names(clinical))) stop("Missing patient ID, PFS or event indicator.")
  clinical <- clinical[!is.na(clinical$FL_ID), fields]
  if (anyDuplicated(clinical$FL_ID) || !all(post$patient_id %in% clinical$FL_ID)) {
    stop("Paired biopsies must match a unique clinical record.")
  }
  clinical <- clinical[match(post$patient_id, clinical$FL_ID), ]
  d <- data.frame(patient_id = post$patient_id, biopsy_timepoint = post$timepoint,
                  delta = delta, group = ifelse(delta > 0, "increase", "decrease"),
                  months = as.numeric(clinical[["PFS (month)"]]), event = as.numeric(clinical$censor))
  if (any(!is.finite(d$months)) || any(d$months < 0) || !all(d$event %in% c(0, 1)) ||
      length(unique(d$group)) != 2L) stop("Invalid PFS times, event indicators or groups.")
  fit <- survival::survfit(survival::Surv(months, event) ~ group, data = d)
  test <- survival::survdiff(survival::Surv(months, event) ~ group, data = d)
  km <- summary(fit, censored = TRUE)
  coords <- data.frame(group = sub("^group=", "", as.character(km$strata)), months = km$time,
                       at_risk = km$n.risk, events = km$n.event, censored = km$n.censor,
                       PFS_probability = km$surv)
  tests <- data.frame(patients = nrow(d), increased = sum(d$group == "increase"),
                      decreased = sum(d$group == "decrease"), chisq = unname(test$chisq),
                      logrank_P = stats::pchisq(test$chisq, 1, lower.tail = FALSE))
  dest <- args[3]
  if (file.exists(dest) && (!dir.exists(dest) || length(list.files(dest, all.files = TRUE, no.. = TRUE)))) {
    stop("Choose a new or empty output directory.")
  }
  if (!dir.exists(dest) && !dir.create(dest, recursive = TRUE)) stop("Cannot create output directory.")
  utils::write.csv(d, file.path(dest, "Fig1C_groups.csv"), row.names = FALSE)
  utils::write.csv(coords, file.path(dest, "Fig1C_KM_coordinates.csv"), row.names = FALSE)
  utils::write.csv(tests, file.path(dest, "Fig1C_logrank.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(dest, "sessionInfo.txt"))
  invisible(tests)
}

if (sys.nframe() == 0L) main()
