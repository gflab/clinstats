test_that("survival_cutoff finds a signal in a strongly separated marker", {
  set.seed(4)
  n <- 400
  marker <- rnorm(n)
  high <- marker > 0
  time <- ifelse(high, stats::runif(n, 0, 0.5), stats::runif(n, 2, 5))
  event <- rep(1, n)
  result <- survival_cutoff(time, event, marker, predict_time = 1)

  expect_s3_class(result, "tbl_df")
  expect_gt(result$youden, 0.8)
  expect_lt(result$cutoff, 1)
  expect_gt(result$tp, 0.8)
  expect_lt(result$fp, 0.1)
})

test_that("survival_cutoff validates its inputs", {
  time <- rexp(50)
  event <- rbinom(50, 1, 0.5)
  marker <- rnorm(50)
  expect_error(survival_cutoff(time, event, marker[-1]), "same length")
  expect_error(
    survival_cutoff(time, rep(2, 50), marker, predict_time = 1),
    "0/1 or logical"
  )
  expect_s3_class(
    survival_cutoff(time, rep(1, 50), marker, predict_time = 1),
    "tbl_df"
  )
})

test_that("survival_response matches a direct Surv call", {
  response <- survival_response(clin_crc, type = "rfs")
  expect_s3_class(response, "Surv")
  expect_equal(
    as.numeric(response[, 1]),
    as.numeric(clin_crc$rfs.delay)
  )
  expect_error(survival_response(clin_crc, type = "none"), "arg")
})

test_that("resample_cox reports selection frequencies and is reproducible", {
  set.seed(5)
  n <- 80
  expression <- matrix(
    rnorm(n * 3),
    nrow = n,
    dimnames = list(NULL, c("g1", "g2", "g3"))
  )
  time <- rexp(n)
  event <- rbinom(n, 1, 0.7)
  first <- resample_cox(expression, time, event, times = 15, seed = 11)
  second <- resample_cox(expression, time, event, times = 15, seed = 11)

  expect_s3_class(first, "tbl_df")
  expect_equal(first, second)
  expect_equal(nrow(first), 3L)
  expect_true(all(first$fraction >= 0 & first$fraction <= 1))
  expect_true(all(first$n_valid == 15))
})
