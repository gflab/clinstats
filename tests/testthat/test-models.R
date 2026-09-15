test_that("cox_table agrees with a direct coxph fit", {
  dat <- toy_survival_data()
  result <- cox_table(dat, time = "time", event = "event", factors = "age")

  reference <- survival::coxph(
    survival::Surv(dat$time, dat$event) ~ dat$age
  )
  expect_equal(
    result$hr_univariable,
    unname(exp(stats::coef(reference)))
  )

  multivariable <- cox_table(
    dat,
    time = "time", event = "event",
    factors = c("age", "stage", "sex")
  )
  expect_equal(nrow(multivariable), 3L)
  expect_false(anyNA(multivariable$hr_multivariable))
})

test_that("cox_table supports screening and a follow-up limit", {
  dat <- toy_survival_data()
  all_variables <- cox_table(
    dat,
    time = "time", event = "event",
    factors = c("age", "stage", "sex")
  )

  screened <- cox_table(
    dat,
    time = "time", event = "event",
    factors = c("age", "stage", "sex"),
    multivariable = "significant"
  )
  significant <- all_variables$univ_p < 0.05
  expect_equal(
    is.na(screened$hr_multivariable),
    !significant
  )
  if (all(significant)) {
    expect_false(anyNA(screened$hr_multivariable))
  }

  limited <- cox_table(
    dat,
    time = "time", event = "event",
    factors = "age", max_time = stats::median(dat$time)
  )
  expect_true(is.finite(limited$hr_univariable))
})

test_that("logistic_table agrees with a direct glm fit", {
  dat <- toy_outcome_data()
  result <- logistic_table(dat, outcome = "outcome", factors = "marker")

  reference <- stats::glm(
    dat$outcome ~ dat$marker,
    family = stats::binomial()
  )
  expect_equal(
    result$or_univariable,
    unname(exp(stats::coef(reference)[2]))
  )
  expect_lt(result$univ_p, 0.01)
})

test_that("logistic_table handles factor levels and validation", {
  dat <- toy_survival_data()
  result <- logistic_table(
    dat,
    outcome = "event",
    factors = c("stage", "sex")
  )
  expect_equal(sort(result$variable), c("sex", "stage"))
  expect_true("III" %in% stats::na.omit(result$level))

  expect_error(logistic_table(dat, outcome = "missing"), "missing")
  expect_error(
    logistic_table(dat, outcome = "time"),
    "positive"
  )
})

test_that("evaluate_model reproduces pROC AUC values", {
  set.seed(3)
  train <- data.frame(
    y = rbinom(100, 1, 0.5),
    x = rnorm(100)
  )
  test <- data.frame(
    y = rbinom(60, 1, 0.5),
    x = rnorm(60)
  )
  result <- evaluate_model(train, test, y ~ x)

  fit <- stats::glm(y ~ x, data = train, family = stats::binomial())
  reference_train <- pROC::roc(
    train$y,
    stats::predict(fit, newdata = train, type = "response"),
    direction = ">",
    quiet = TRUE
  )
  reference_test <- pROC::roc(
    test$y,
    stats::predict(fit, newdata = test, type = "response"),
    direction = ">",
    quiet = TRUE
  )
  expect_equal(result$auc_train, as.numeric(pROC::auc(reference_train)))
  expect_equal(result$auc_test, as.numeric(pROC::auc(reference_test)))
  expect_true(result$delong_p >= 0 && result$delong_p <= 1)
})

test_that("column validation reports the offending argument", {
  dat <- toy_survival_data()
  expect_error(
    cox_table(dat, time = "nope", event = "event"),
    "time"
  )
  expect_error(
    cox_table(dat, time = "event", event = "time", factors = "missing"),
    "missing"
  )
})
