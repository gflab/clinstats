test_that("roc_summary matches pROC reference values", {
  dat <- toy_outcome_data()
  reference <- pROC::roc(dat$outcome, dat$marker, quiet = TRUE)
  result <- roc_summary(dat$marker, dat$outcome)

  expect_s3_class(result, "tbl_df")
  expect_equal(result$auc, as.numeric(pROC::auc(reference)))
  expect_equal(result$n, nrow(dat))
  expect_equal(result$n_event, sum(dat$outcome == 1))
  expect_true(result$sensitivity >= 0 && result$sensitivity <= 1)
  expect_true(result$specificity >= 0 && result$specificity <= 1)
})

test_that("roc_summary handles logical and labelled outcomes", {
  dat <- toy_outcome_data()
  logical_result <- roc_summary(dat$marker, dat$outcome == 1)
  expect_equal(logical_result$auc, roc_summary(dat$marker, dat$outcome)$auc)

  labelled <- ifelse(dat$outcome == 1, "case", "control")
  expect_error(roc_summary(dat$marker, labelled), "positive")
  labelled_result <- roc_summary(dat$marker, labelled, positive = "case")
  expect_equal(labelled_result$auc, logical_result$auc)
})

test_that("roc_summary validates its inputs", {
  dat <- toy_outcome_data()
  expect_error(roc_summary(dat$marker, dat$outcome[-1]), "same length")
  expect_error(roc_summary("a", dat$outcome), "numeric")
})

test_that("odds_ratio reproduces a hand-computed odds ratio", {
  dat <- toy_outcome_data()
  result <- odds_ratio(dat$marker, dat$outcome)

  group <- dat$marker > result$cutoff
  table_2x2 <- table(group, dat$outcome)
  manual <- (table_2x2["TRUE", "1"] * table_2x2["FALSE", "0"]) /
    (table_2x2["TRUE", "0"] * table_2x2["FALSE", "1"])
  expect_equal(result$odds_ratio, unname(manual))
  expect_true(result$ci_lower < result$odds_ratio)
  expect_true(result$ci_upper > result$odds_ratio)
})

test_that("odds_ratio accepts two-level and factor inputs", {
  dat <- toy_outcome_data()
  binary_marker <- dat$marker > stats::median(dat$marker)
  expect_equal(
    odds_ratio(binary_marker, dat$outcome)$odds_ratio,
    odds_ratio(binary_marker, dat$outcome, cutoff = 0.5)$odds_ratio
  )
  grouped <- factor(ifelse(dat$marker > 0.5, "high", "low"), levels = c("low", "high"))
  result <- odds_ratio(grouped, dat$outcome)
  expect_equal(result$cutoff, "high")
  expect_true(result$odds_ratio > 1)
})

test_that("odds_ratio rejects thresholds that do not split the data", {
  dat <- toy_outcome_data()
  expect_error(
    odds_ratio(dat$marker, dat$outcome, cutoff = 100),
    "two groups"
  )
})
