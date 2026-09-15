# Internal helpers --------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

check_numeric <- function(x, arg) {
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg {arg}} must be a numeric vector, not {.cls {class(x)[1]}}.")
  }
  invisible(x)
}

check_length <- function(x, y, x_arg, y_arg) {
  if (length(x) != length(y)) {
    cli::cli_abort(c(
      "{.arg {x_arg}} and {.arg {y_arg}} must have the same length.",
      x = "{.arg {x_arg}} has length {length(x)}.",
      x = "{.arg {y_arg}} has length {length(y)}."
    ))
  }
  invisible(x)
}

check_columns <- function(data, columns, arg = "data") {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg {arg}} must be a data frame, not {.cls {class(data)[1]}}.")
  }
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "{.arg {arg}} is missing the required column{?s} {.val {missing}}.",
      i = "Available columns: {.val {names(data)}}."
    ))
  }
  invisible(data)
}

# Resolve a binary outcome to a 0/1 integer vector. Logical outcomes are
# coded TRUE = 1. Numeric 0/1 outcomes use 1 as the event. Anything else
# requires an explicit `positive` value so the event class is never guessed.
resolve_binary_outcome <- function(outcome, positive = NULL, arg = "outcome") {
  if (is.logical(outcome)) {
    positive <- positive %||% TRUE
  } else {
    values <- unique(outcome[!is.na(outcome)])
    if (is.numeric(outcome) && all(values %in% c(0, 1))) {
      positive <- positive %||% 1
    } else if (is.null(positive)) {
      cli::cli_abort(c(
        "{.arg {arg}} is not a 0/1 or logical vector.",
        i = "Supply {.arg positive} to name the event level explicitly."
      ))
    }
  }
  values <- unique(outcome[!is.na(outcome)])
  if (!positive %in% values) {
    cli::cli_abort("{.val {positive}} is not a value of {.arg {arg}}.")
  }
  if (length(values) != 2) {
    cli::cli_abort("{.arg {arg}} must have exactly two distinct non-missing values.")
  }
  as.integer(outcome == positive)
}

# Wald confidence interval on the link scale.
wald_ci <- function(estimate, se, level) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  rbind(lower = estimate - z * se, upper = estimate + z * se)
}

# Resolve an event indicator to 0/1 without requiring both values to occur:
# a sample in which every observation has an event is a valid input to the
# Kaplan-Meier estimator.
resolve_event_indicator <- function(event, arg = "event") {
  if (is.logical(event)) {
    return(as.integer(event))
  }
  if (is.numeric(event) && all(event[!is.na(event)] %in% c(0, 1))) {
    return(as.integer(event))
  }
  cli::cli_abort("{.arg {arg}} must be a 0/1 or logical vector.")
}

# Tidy the coefficient table of a fitted Cox or logistic model. Returns one
# row per model coefficient with Wald confidence limits on the exponentiated
# (hazard ratio / odds ratio) scale.
tidy_regression <- function(fit, ci_level = 0.95) {
  co <- summary(fit)$coefficients
  co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  beta <- co[, grep("^(coef|Estimate)$", colnames(co))[1]]
  se <- co[, grep("^(se\\(coef\\)|Std\\. Error)$", colnames(co))[1]]
  p_col <- grep("^Pr\\(", colnames(co), value = TRUE)[1]
  ci <- wald_ci(beta, se, ci_level)
  tibble::tibble(
    coefficient = rownames(co),
    estimate = exp(beta),
    ci_lower = exp(ci["lower", ]),
    ci_upper = exp(ci["upper", ]),
    p_value = co[, p_col]
  )
}

# Build the working data frame used by the regression tables. Variables are
# copied under short synthetic names (`v01`, `v02`, ...) so that model
# coefficients can be mapped back to variables and factor levels without
# name-prefix ambiguity.
regression_work <- function(data, factors, response = NULL) {
  work <- data.frame(row.names = seq_len(nrow(data)))
  if (!is.null(response)) {
    work$.response <- response
  }
  for (i in seq_along(factors)) {
    work[[sprintf("v%02d", i)]] <- data[[factors[i]]]
  }
  work
}

# Map the coefficient rows of a model fitted on `regression_work()` output
# back to variable names, factor levels, and readable terms.
map_coefficients <- function(tidy, factors) {
  index <- as.integer(substr(tidy$coefficient, 2L, 3L))
  tidy$variable <- factors[index]
  tidy$level <- substring(tidy$coefficient, 4L)
  tidy$level[tidy$level == ""] <- NA_character_
  tidy$term <- ifelse(
    is.na(tidy$level),
    tidy$variable,
    paste0(tidy$variable, tidy$level)
  )
  tidy
}

# Variable-level p-value used to select terms for the multivariable model.
variable_p_value <- function(fit, type) {
  table <- if (type == "cox") {
    stats::anova(fit)
  } else {
    stats::anova(fit, test = "LRT")
  }
  p_col <- grep("^Pr\\(", colnames(table), value = TRUE)[1]
  table[nrow(table), p_col]
}
