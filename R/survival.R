#' Survival cut point optimising the Youden index
#'
#' Estimates Kaplan-Meier survival above and below each candidate cut point
#' of a marker at `predict_time` and returns the cut point that maximises the
#' Youden index (sensitivity + specificity - 1), the same criterion that the
#' long-unmaintained `survivalROC` package uses, computed here directly from
#' [survival::survfit()]. Candidates are drawn from a quantile grid to keep
#' the computation bounded.
#'
#' @param time Numeric follow-up times.
#' @param event Event indicator (0/1, logical, or a two-level value).
#' @param marker Numeric marker to split.
#' @param predict_time Time point at which sensitivity and specificity are
#'   evaluated.
#' @param grid_size Number of candidate cut points evaluated between the 1st
#'   and 99th percentile of the marker.
#'
#' @return A tibble with one row: the optimal `cutoff`, the Youden index
#'   (`youden`), the time-dependent true positive rate (`tp`) and false
#'   positive rate (`fp`) at that cut point, and the number of observations.
#'
#' @export
#' @examples
#' set.seed(1)
#' time <- rexp(100)
#' event <- rbinom(100, 1, 0.6)
#' marker <- rnorm(100)
#' survival_cutoff(time, event, marker, predict_time = 1)
survival_cutoff <- function(time, event, marker, predict_time,
                            grid_size = 200) {
  check_numeric(time, "time")
  check_numeric(marker, "marker")
  check_length(time, marker, "time", "marker")
  check_length(time, event, "time", "event")
  if (!is.numeric(predict_time) || length(predict_time) != 1) {
    cli::cli_abort("{.arg predict_time} must be a single number.")
  }
  ok <- !is.na(time) & !is.na(event) & !is.na(marker)
  time <- time[ok]
  marker <- marker[ok]
  event <- resolve_event_indicator(event[ok], arg = "event")
  if (length(time) < 5) {
    cli::cli_abort("At least five complete observations are required.")
  }

  probs <- seq(0.01, 0.99, length.out = grid_size)
  candidates <- unique(stats::quantile(marker, probs = probs, na.rm = TRUE))
  km_surv <- function(index) {
    fit <- survival::survfit(survival::Surv(time[index], event[index]) ~ 1)
    s <- summary(fit, times = predict_time, extend = TRUE)$surv
    if (length(s) == 0) NA_real_ else s[1]
  }

  evaluate <- vapply(candidates, function(cut) {
    high <- marker > cut
    if (all(high) || !any(high)) {
      return(c(tp = NA_real_, fp = NA_real_))
    }
    c(tp = 1 - km_surv(high), fp = 1 - km_surv(!high))
  }, numeric(2))

  youden <- evaluate["tp", ] - evaluate["fp", ]
  best <- which.max(youden)
  if (length(best) == 0 || is.na(youden[best])) {
    cli::cli_abort("No cut point could be evaluated at {.arg predict_time}.")
  }
  tibble::tibble(
    cutoff = candidates[best],
    youden = youden[best],
    tp = evaluate["tp", best],
    fp = evaluate["fp", best],
    predict_time = predict_time,
    n = length(time),
    n_event = sum(event == 1)
  )
}

#' Build a survival response from a clinical data frame
#'
#' @param data Data frame such as the included [clin_crc] dataset.
#' @param type Endpoint prefix: `"rfs"`, `"os"`, or `"dfs"`. The columns
#'   `<type>.delay` and `<type>.event` must exist.
#'
#' @return A [survival::Surv()] object.
#'
#' @export
#' @examples
#' survival_response(clin_crc, type = "rfs")
survival_response <- function(data, type = c("rfs", "os", "dfs")) {
  type <- match.arg(type)
  time_col <- paste0(type, ".delay")
  event_col <- paste0(type, ".event")
  check_columns(data, c(time_col, event_col))
  survival::Surv(data[[time_col]], data[[event_col]] == 1)
}

#' Resampling-based Cox screening of candidate markers
#'
#' Repeatedly draws subsamples of the samples, fits a univariable Cox model
#' for every candidate column, and reports how often each candidate reached
#' `p_threshold`. This is the resampling screening previously provided by
#' `calc_resamp_cox()`, reported as per-marker selection frequencies instead
#' of a raw p-value matrix.
#'
#' The function uses the `future` framework. Set a parallel plan before
#' calling it to use several workers, for example
#' `future::plan(future::multisession)`.
#'
#' @param expression Numeric matrix or data frame with samples in rows and
#'   candidate markers in columns.
#' @param time,event Follow-up time and event indicator, one value per sample.
#' @param ratio Fraction of samples drawn in each resample.
#' @param times Number of resamples.
#' @param p_threshold Significance threshold that counts as a selection.
#' @param seed Optional integer seed for reproducible resampling.
#'
#' @return A tibble with one row per candidate marker: the number of valid
#'   resamples, the number of resamples in which the marker was significant,
#'   the selection fraction, and the median p-value.
#'
#' @export
#' @examples
#' set.seed(1)
#' expr <- matrix(rnorm(60 * 5), nrow = 60, dimnames = list(NULL, paste0("g", 1:5)))
#' time <- rexp(60)
#' event <- rbinom(60, 1, 0.7)
#' resample_cox(expr, time, event, times = 20, seed = 1)
resample_cox <- function(expression, time, event, ratio = 0.8, times = 1000,
                         p_threshold = 0.05, seed = NULL) {
  expression <- as.matrix(expression)
  if (!is.numeric(expression)) {
    cli::cli_abort("{.arg expression} must be numeric.")
  }
  check_length(time, event, "time", "event")
  if (nrow(expression) != length(time)) {
    cli::cli_abort(c(
      "{.arg expression} must have one row per sample.",
      x = "{.arg expression} has {nrow(expression)} rows.",
      x = "{.arg time} has length {length(time)}."
    ))
  }
  if (is.null(colnames(expression))) {
    cli::cli_abort("{.arg expression} must have column names for the markers.")
  }
  if (ratio <= 0 || ratio > 1 || times < 1) {
    cli::cli_abort("{.arg ratio} must be in (0, 1] and {.arg times} >= 1.")
  }
  event <- resolve_event_indicator(event, arg = "event")
  ok <- !is.na(time) & !is.na(event)
  expression <- expression[ok, , drop = FALSE]
  time <- time[ok]
  event <- event[ok]
  n <- nrow(expression)
  size <- max(2, floor(n * ratio))

  screen <- function(i) {
    index <- sample.int(n, size)
    apply(expression[index, , drop = FALSE], 2, function(g) {
      tryCatch(
        {
          fit <- survival::coxph(survival::Surv(time[index], event[index]) ~ g)
          summary(fit)$coefficients[, "Pr(>|z|)"]
        },
        error = function(e) NA_real_
      )
    })
  }

  run <- function() {
    future.apply::future_lapply(seq_len(times), screen, future.seed = TRUE)
  }
  draws <- if (is.null(seed)) run() else withr::with_seed(seed, run())
  p_values <- do.call(cbind, draws)
  n_valid <- rowSums(!is.na(p_values))

  tibble::tibble(
    gene = colnames(expression),
    n_valid = n_valid,
    n_significant = rowSums(p_values < p_threshold, na.rm = TRUE),
    fraction = rowSums(p_values < p_threshold, na.rm = TRUE) / n_valid,
    median_p = apply(p_values, 1, stats::median, na.rm = TRUE)
  )
}
