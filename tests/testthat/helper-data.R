set.seed(2026)

toy_outcome_data <- function(n = 120, seed = 1) {
  set.seed(seed)
  marker <- c(rnorm(n / 2, 0), rnorm(n / 2, 1.2))
  tibble::tibble(
    marker = marker,
    outcome = rep(c(0L, 1L), each = n / 2)
  )
}

toy_survival_data <- function(n = 120, seed = 2) {
  set.seed(seed)
  age <- rnorm(n, 60, 8)
  stage <- factor(sample(c("II", "III"), n, replace = TRUE))
  time <- rexp(n, rate = 0.08 * exp(stage == "III") * exp(age - 60) / 20)
  tibble::tibble(
    time = time,
    event = rbinom(n, 1, 0.7),
    age = age,
    stage = stage,
    sex = factor(sample(c("F", "M"), n, replace = TRUE))
  )
}
