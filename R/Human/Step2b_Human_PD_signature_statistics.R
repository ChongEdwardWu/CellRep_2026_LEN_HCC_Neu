#!/usr/bin/env Rscript

# The bootstrap samples patients independently with replacement within each
# response group. For an odd sample size the median is the middle order statistic.
bootstrap_median <- function(x) {
  n <- length(x)
  if (n %% 2L == 0L) stop("The median distribution here requires an odd group size.")
  values <- sort(unique(x))
  empirical_cdf <- cumsum(tabulate(match(x, values), nbins = length(values))) / n
  # Median <= v when at least (n + 1)/2 draws are <= v.
  cdf <- stats::pbinom((n - 1L) / 2L, n, empirical_cdf, lower.tail = FALSE)
  list(values = values, mass = diff(c(0, cdf)))
}

median_difference_ci <- function(pd, non_pd) {
  a <- bootstrap_median(pd)
  b <- bootstrap_median(non_pd)
  differences <- as.vector(outer(a$values, b$values, "-"))
  weights <- as.vector(outer(a$mass, b$mass, "*"))
  index <- order(differences)
  cdf <- cumsum(weights[index])
  vapply(c(0.025, 0.975), function(p) differences[index[which(cdf >= p)[1L]]], numeric(1))
}

signature_statistics <- function(signature, data) {
  rows <- data[data$signature == signature, ]
  pd <- rows$delta[rows$group == "PD"]
  non_pd <- rows$delta[rows$group == "non-PD"]
  p <- if (length(unique(c(pd, non_pd))) == 1L) 1 else
    stats::wilcox.test(pd, non_pd, alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value
  # Pairwise wins minus losses per PD/non-PD pair. A tied pair contributes zero.
  r <- mean(sign(outer(pd, non_pd, "-")))
  ci <- median_difference_ci(pd, non_pd)
  data.frame(signature = signature, n_PD = length(pd), n_nonPD = length(non_pd),
             mean_PD = mean(pd), mean_nonPD = mean(non_pd), diff_mean = mean(pd) - mean(non_pd),
             median_PD = median(pd), median_nonPD = median(non_pd),
             diff_median = median(pd) - median(non_pd),
             prop_up_PD = mean(pd > 0), prop_up_nonPD = mean(non_pd > 0),
             p_wilcox = p, r_rank_biserial = r, ci_med_low = ci[1], ci_med_high = ci[2])
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Usage: Rscript Step2b_Human_PD_signature_statistics.R",
                 "<PD_delta_long_by_patient_signature.csv> <output_dir>")
  if (length(args) == 1L && args %in% c("-h", "--help")) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  if (length(args) != 2L) stop(usage, call. = FALSE)
  data <- utils::read.csv(args[1], stringsAsFactors = FALSE)
  required <- c("patient_id", "group", "signature", "delta")
  if (!all(required %in% names(data)) || anyNA(data[required]) ||
      !is.numeric(data$delta) || any(!is.finite(data$delta))) {
    stop("Expected Step 2 patient IDs, groups, signatures and finite numeric paired changes.")
  }
  groups <- unique(data[c("patient_id", "group")])
  signatures <- unique(data$signature)
  if (any(!nzchar(data$patient_id)) || any(!nzchar(signatures)) || anyDuplicated(groups$patient_id) ||
      !setequal(groups$group, c("PD", "non-PD")) || sum(groups$group == "PD") != 3L ||
      sum(groups$group == "non-PD") != 9L || length(signatures) != 10L ||
      any(table(data$patient_id, data$signature) != 1L)) {
    stop("Figure S1B requires one change per patient and signature for three PD and nine non-PD patients.")
  }
  stats <- do.call(rbind, lapply(signatures, signature_statistics, data = data))
  stats$q_fdr <- stats::p.adjust(stats$p_wilcox, method = "BH")
  stats <- stats[order(-abs(stats$diff_median), stats$signature), ]
  out <- args[2]
  if (file.exists(out) && (!dir.exists(out) || length(list.files(out, all.files = TRUE, no.. = TRUE)))) {
    stop("Output directory must be new or empty.")
  }
  if (!dir.exists(out) && !dir.create(out, recursive = TRUE)) stop("Could not create output directory.")
  utils::write.csv(stats, file.path(out, "S1B_STable_PD_delta_stats_ALLsignatures.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  invisible(stats)
}

if (sys.nframe() == 0L) main()
