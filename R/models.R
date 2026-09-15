#' Regression tables for clinical variables
#'
#' `cox_table()` and `logistic_table()` fit each variable on its own, fit one
#' multivariable model, and return a single tidy table with the univariable
#' and multivariable estimates side by side. Numeric variables contribute one
#' row with a per-unit estimate; factors contribute one row per level against
#' the reference level.
#'
#' The multivariable model includes all requested variables by default. With
#' `multivariable = "significant"` only variables whose univariable
#' likelihood-ratio test is below `significance` are included, which
#' reproduces the screening rule of the original helper functions.
#'
#' @param data Data frame holding the variables.
#' @param time,event Column names of the follow-up time and event indicator
#'   used by `cox_table()`.
#' @param outcome Column name of the binary outcome used by
#'   `logistic_table()`.
#' @param factors Character vector of variable names to analyse. Defaults to
#'   every column of `data` except the time and event columns.
#' @param max_time Optional follow-up limit for `cox_table()`. Times beyond
#'   the limit are censored at the limit.
#' @param multivariable Which variables enter the multivariable model:
#'   `"all"` (default), `"significant"`, or `"none"`.
#' @param significance Significance level used by
#'   `multivariable = "significant"`.
#' @param ci_level Confidence level of the reported intervals.
#'
#' @return A tibble with one row per term and the columns `term`, `variable`,
#'   `level`, `n`, `n_event`, the univariable estimate (`hr_univariable` or
#'   `or_univariable`) with interval and p-value, and the corresponding
#'   multivariable columns.
#'
#' @export
#' @examples
#' set.seed(42)
#' dat <- data.frame(
#'   time = rexp(60, 0.1),
#'   event = rbinom(60, 1, 0.7),
#'   age = rnorm(60, 60, 8),
#'   stage = factor(sample(c("II", "III"), 60, replace = TRUE))
#' )
#' cox_table(dat, time = "time", event = "event")
cox_table <- function(data, time, event, factors = NULL, max_time = NULL,
                      multivariable = c("all", "significant", "none"),
                      significance = 0.05, ci_level = 0.95) {
  multivariable <- match.arg(multivariable)
  check_columns(data, c(time, event))
  if (is.null(factors)) factors <- setdiff(names(data), c(time, event))
  check_columns(data, factors)
  if (length(factors) < 1) {
    cli::cli_abort("At least one variable is required in {.arg factors}.")
  }

  tvec <- data[[time]]
  evec <- data[[event]]
  check_numeric(tvec, "time")
  if (!is.null(max_time)) {
    evec <- ifelse(tvec > max_time, 0, evec)
    tvec <- pmin(tvec, max_time)
  }
  work <- regression_work(data, factors)
  work$.time <- tvec
  work$.event <- evec

  univ_fits <- lapply(seq_along(factors), function(i) {
    survival::coxph(
      stats::as.formula(paste(
        "survival::Surv(.time, .event) ~", sprintf("v%02d", i)
      )),
      data = work
    )
  })
  univ <- do.call(rbind, lapply(univ_fits, function(fit) {
    tidy_regression(fit, ci_level = ci_level)
  }))
  univ <- map_coefficients(univ, factors)
  univ$n <- vapply(seq_along(factors), function(i) {
    sum(!is.na(data[[factors[i]]]) & !is.na(tvec) & !is.na(evec))
  }, integer(1))[match(univ$variable, factors)]
  univ$n_event <- vapply(seq_along(factors), function(i) {
    ok <- !is.na(data[[factors[i]]]) & !is.na(tvec) & !is.na(evec)
    sum(evec[ok] == 1)
  }, integer(1))[match(univ$variable, factors)]

  multi <- NULL
  if (multivariable != "none") {
    keep <- rep(TRUE, length(factors))
    if (multivariable == "significant") {
      p <- vapply(univ_fits, variable_p_value, numeric(1), type = "cox")
      keep <- p < significance
    }
    if (any(keep)) {
      multi_fit <- survival::coxph(
        stats::as.formula(paste(
          "survival::Surv(.time, .event) ~",
          paste(sprintf("v%02d", which(keep)), collapse = " + ")
        )),
        data = work
      )
      multi <- map_coefficients(tidy_regression(multi_fit, ci_level), factors)
    }
  }

  index <- if (is.null(multi)) rep(NA_integer_, nrow(univ)) else {
    match(univ$term, multi$term)
  }
  tibble::tibble(
    term = univ$term,
    variable = univ$variable,
    level = univ$level,
    n = univ$n,
    n_event = univ$n_event,
    hr_univariable = univ$estimate,
    univ_ci_lower = univ$ci_lower,
    univ_ci_upper = univ$ci_upper,
    univ_p = univ$p_value,
    hr_multivariable = if (is.null(multi)) NA_real_ else multi$estimate[index],
    multi_ci_lower = if (is.null(multi)) NA_real_ else multi$ci_lower[index],
    multi_ci_upper = if (is.null(multi)) NA_real_ else multi$ci_upper[index],
    multi_p = if (is.null(multi)) NA_real_ else multi$p_value[index]
  )
}

#' @rdname cox_table
#' @param positive Value of `outcome` that represents the event.
#' @export
#' @examples
#' set.seed(42)
#' dat <- data.frame(
#'   response = rbinom(60, 1, 0.4),
#'   age = rnorm(60, 60, 8),
#'   sex = factor(sample(c("F", "M"), 60, replace = TRUE))
#' )
#' logistic_table(dat, outcome = "response")
logistic_table <- function(data, outcome, factors = NULL,
                           positive = NULL,
                           multivariable = c("all", "significant", "none"),
                           significance = 0.05, ci_level = 0.95) {
  multivariable <- match.arg(multivariable)
  check_columns(data, outcome)
  if (is.null(factors)) factors <- setdiff(names(data), outcome)
  check_columns(data, factors)
  if (length(factors) < 1) {
    cli::cli_abort("At least one variable is required in {.arg factors}.")
  }

  y <- resolve_binary_outcome(data[[outcome]], positive, arg = "outcome")
  work <- regression_work(data, factors, y)

  univ_fits <- lapply(seq_along(factors), function(i) {
    stats::glm(
      stats::as.formula(paste(".response ~", sprintf("v%02d", i))),
      family = stats::binomial(),
      data = work
    )
  })
  univ <- do.call(rbind, lapply(univ_fits, function(fit) {
    tidy_regression(fit, ci_level = ci_level)
  }))
  univ <- map_coefficients(univ, factors)
  univ$n <- vapply(seq_along(factors), function(i) {
    sum(!is.na(data[[factors[i]]]) & !is.na(data[[outcome]]))
  }, integer(1))[match(univ$variable, factors)]
  univ$n_event <- vapply(seq_along(factors), function(i) {
    ok <- !is.na(data[[factors[i]]]) & !is.na(data[[outcome]])
    sum(y[ok] == 1)
  }, integer(1))[match(univ$variable, factors)]

  multi <- NULL
  if (multivariable != "none") {
    keep <- rep(TRUE, length(factors))
    if (multivariable == "significant") {
      p <- vapply(univ_fits, variable_p_value, numeric(1), type = "glm")
      keep <- p < significance
    }
    if (any(keep)) {
      multi_fit <- stats::glm(
        stats::as.formula(paste(
          ".response ~",
          paste(sprintf("v%02d", which(keep)), collapse = " + ")
        )),
        family = stats::binomial(),
        data = work
      )
      multi <- map_coefficients(tidy_regression(multi_fit, ci_level), factors)
    }
  }

  index <- if (is.null(multi)) rep(NA_integer_, nrow(univ)) else {
    match(univ$term, multi$term)
  }
  tibble::tibble(
    term = univ$term,
    variable = univ$variable,
    level = univ$level,
    n = univ$n,
    n_event = univ$n_event,
    or_univariable = univ$estimate,
    univ_ci_lower = univ$ci_lower,
    univ_ci_upper = univ$ci_upper,
    univ_p = univ$p_value,
    or_multivariable = if (is.null(multi)) NA_real_ else multi$estimate[index],
    multi_ci_lower = if (is.null(multi)) NA_real_ else multi$ci_lower[index],
    multi_ci_upper = if (is.null(multi)) NA_real_ else multi$ci_upper[index],
    multi_p = if (is.null(multi)) NA_real_ else multi$p_value[index]
  )
}

#' Evaluate a prediction model on training and test data
#'
#' Fits a logistic regression model on `train`, evaluates discrimination on
#' both `train` and `test`, and compares the two areas under the curve with
#' DeLong's test for unpaired samples.
#'
#' @param train,test Data frames with the model variables.
#' @param formula Model formula, for example `outcome ~ marker1 + marker2`.
#' @param positive Value of the outcome that represents the event.
#' @param ci_level Confidence level of the reported AUC intervals.
#'
#' @return A tibble with one row: sample sizes and event counts of both data
#'   sets, the training and test AUC with confidence intervals, and the
#'   DeLong p-value comparing them.
#'
#' @export
#' @examples
#' set.seed(1)
#' train <- data.frame(y = rbinom(80, 1, 0.5), x = rnorm(80))
#' test <- data.frame(y = rbinom(40, 1, 0.5), x = rnorm(40))
#' evaluate_model(train, test, y ~ x)
evaluate_model <- function(train, test, formula, positive = NULL,
                           ci_level = 0.95) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort("{.arg formula} must be a model formula.")
  }
  outcome <- all.vars(formula[[2]])
  if (length(outcome) != 1) {
    cli::cli_abort("The left-hand side of {.arg formula} must be a single variable.")
  }
  check_columns(train, c(outcome, all.vars(formula[[3]])), arg = "train")
  check_columns(test, c(outcome, all.vars(formula[[3]])), arg = "test")

  fit <- stats::glm(formula, data = train, family = stats::binomial())
  y_train <- resolve_binary_outcome(train[[outcome]], positive)
  y_test <- resolve_binary_outcome(test[[outcome]], positive)
  p_train <- stats::predict(fit, newdata = train, type = "response")
  p_test <- stats::predict(fit, newdata = test, type = "response")

  # Predicted probabilities must keep a fixed direction: larger values mean
  # a higher event probability, so an AUC below 0.5 honestly reports a model
  # that ranks the test data in the wrong direction.
  roc_train <- pROC::roc(y_train, p_train, direction = ">", quiet = TRUE)
  roc_test <- pROC::roc(y_test, p_test, direction = ">", quiet = TRUE)
  ci_train <- as.numeric(pROC::ci.auc(roc_train, conf.level = ci_level))
  ci_test <- as.numeric(pROC::ci.auc(roc_test, conf.level = ci_level))
  delong <- tryCatch(
    stats::p.adjust(
      pROC::roc.test(roc_train, roc_test, method = "delong", paired = FALSE)$p.value,
      method = "none"
    ),
    error = function(e) NA_real_
  )

  tibble::tibble(
    n_train = length(y_train),
    n_event_train = sum(y_train == 1),
    auc_train = as.numeric(pROC::auc(roc_train)),
    train_ci_lower = ci_train[1],
    train_ci_upper = ci_train[3],
    n_test = length(y_test),
    n_event_test = sum(y_test == 1),
    auc_test = as.numeric(pROC::auc(roc_test)),
    test_ci_lower = ci_test[1],
    test_ci_upper = ci_test[3],
    delong_p = as.numeric(delong)
  )
}
