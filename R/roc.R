#' Summarise the diagnostic performance of a numeric marker
#'
#' Computes the area under the receiver operating characteristic curve (AUC)
#' with a confidence interval and the optimal operating point of a numeric
#' marker against a binary outcome. A tibble with one row is returned, so the
#' result can be combined across markers or subgroups.
#'
#' @param marker Numeric vector of marker values.
#' @param outcome Binary outcome. Logical vectors are coded `TRUE` as the
#'   event; numeric 0/1 vectors use 1 as the event. Any other coding requires
#'   `positive` to name the event level.
#' @param positive Value of `outcome` that represents the event.
#' @param direction Direction of the comparison passed on to
#'   [pROC::roc()]. `"auto"` lets \pkg{pROC} choose the direction from the
#'   data; `"<"` and `">"` fix it explicitly.
#' @param best_method Method used to choose the optimal operating point,
#'   `"youden"` (default) or `"closest.topleft"`.
#' @param ci_level Confidence level for the AUC interval.
#' @param power If `TRUE`, add the power of the AUC test as reported by
#'   [pROC::power.roc.test()]. The value may be `NA` when power cannot be
#'   estimated for the sample.
#' @param alternative Alternative hypothesis used for the power estimate.
#'
#' @return A tibble with one row: sample size, number of events, AUC and its
#'   confidence interval, the optimal threshold, and the sensitivity,
#'   specificity, positive and negative predictive values, and accuracy at
#'   that threshold.
#'
#' @seealso [odds_ratio()] for the odds ratio at the optimal threshold.
#' @export
#' @examples
#' set.seed(1)
#' marker <- c(rnorm(40, 0), rnorm(40, 1))
#' outcome <- rep(c(0, 1), each = 40)
#' roc_summary(marker, outcome)
roc_summary <- function(marker, outcome, positive = NULL,
                        direction = "auto",
                        best_method = c("youden", "closest.topleft"),
                        ci_level = 0.95,
                        power = FALSE,
                        alternative = c("one.sided", "two.sided")) {
  best_method <- match.arg(best_method)
  alternative <- match.arg(alternative)
  check_numeric(marker, "marker")
  y <- resolve_binary_outcome(outcome, positive)
  check_length(marker, y, "marker", "outcome")
  ok <- !is.na(marker) & !is.na(y)
  if (sum(ok) < 2 || length(unique(y[ok])) < 2) {
    cli::cli_abort("At least one observation from each outcome class is required.")
  }

  roc <- pROC::roc(y[ok], marker[ok], direction = direction, quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc, conf.level = ci_level))
  point <- pROC::coords(
    roc, "best",
    best.method = best_method,
    ret = c(
      "threshold", "sensitivity", "specificity",
      "ppv", "npv", "accuracy"
    ),
    transpose = FALSE
  )
  point <- point[1, , drop = FALSE]

  out <- tibble::tibble(
    n = sum(ok),
    n_event = sum(y[ok] == 1),
    auc = as.numeric(pROC::auc(roc)),
    auc_ci_lower = auc_ci[1],
    auc_ci_upper = auc_ci[3],
    direction = roc$direction,
    threshold = point$threshold,
    sensitivity = point$sensitivity,
    specificity = point$specificity,
    ppv = point$ppv,
    npv = point$npv,
    accuracy = point$accuracy
  )

  if (isTRUE(power)) {
    out$power <- tryCatch(
      pROC::power.roc.test(roc, alternative = alternative)$power,
      error = function(e) NA_real_
    )
  }
  out
}

#' Odds ratio for a marker with a binary outcome
#'
#' Binarises a marker and reports the odds ratio for the event in the high
#' group relative to the low group. Continuous markers are split at the
#' Youden-optimal ROC threshold unless `cutoff` is supplied; two-level
#' markers and factors are used as they are.
#'
#' @inheritParams roc_summary
#' @param x Marker or factor to binarise. Continuous numeric values are split
#'   at the optimal threshold; numeric, logical, character, and factor inputs
#'   with two distinct values are used directly.
#' @param cutoff Optional threshold supplied by the user. Values greater than
#'   `cutoff` form the high group.
#'
#' @return A tibble with one row: the threshold that was used, sample size,
#'   number of events, group sizes, the odds ratio with a Wald confidence
#'   interval, and the Wald p-value.
#'
#' @export
#' @examples
#' set.seed(1)
#' marker <- c(rnorm(40, 0), rnorm(40, 1))
#' outcome <- rep(c(0, 1), each = 40)
#' odds_ratio(marker, outcome)
odds_ratio <- function(x, outcome, positive = NULL, cutoff = NULL,
                       ci_level = 0.95) {
  if (is.logical(x)) {
    x <- as.integer(x)
  }
  y <- resolve_binary_outcome(outcome, positive)
  check_length(x, y, "x", "outcome")
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  if (length(y) < 4 || length(unique(y)) < 2) {
    cli::cli_abort("At least one observation from each outcome class is required.")
  }

  cutoff_used <- cutoff
  if (is.null(cutoff) && is.numeric(x) && length(unique(x)) > 2) {
    roc <- pROC::roc(y, x, quiet = TRUE)
    cutoff_used <- as.numeric(pROC::coords(
      roc, "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    )$threshold)[1]
    group <- x > cutoff_used
  } else if (is.null(cutoff)) {
    values <- if (is.factor(x)) {
      levels(x)
    } else {
      sort(unique(as.character(x)))
    }
    if (length(values) != 2) {
      cli::cli_abort(c(
        "{.arg x} must have exactly two distinct values for a direct odds ratio.",
        i = "Provide {.arg cutoff} or a continuous numeric marker."
      ))
    }
    group <- as.character(x) == values[2]
    cutoff_used <- values[2]
  } else {
    check_numeric(x, "x")
    group <- x > cutoff
  }

  if (length(unique(group)) < 2) {
    cli::cli_abort("The threshold does not separate the data into two groups.")
  }

  fit <- stats::glm(y ~ group, family = stats::binomial())
  co <- summary(fit)$coefficients
  ci <- wald_ci(co[2, "Estimate"], co[2, "Std. Error"], ci_level)

  tibble::tibble(
    cutoff = cutoff_used,
    n = length(y),
    n_event = sum(y == 1),
    n_high = sum(group),
    n_low = sum(!group),
    odds_ratio = exp(co[2, "Estimate"]),
    ci_lower = exp(ci["lower", 1]),
    ci_upper = exp(ci["upper", 1]),
    p_value = co[2, "Pr(>|z|)"]
  )
}
