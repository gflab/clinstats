# Regression tables for clinical variables

`cox_table()` and `logistic_table()` fit each variable on its own, fit
one multivariable model, and return a single tidy table with the
univariable and multivariable estimates side by side. Numeric variables
contribute one row with a per-unit estimate; factors contribute one row
per level against the reference level.

## Usage

``` r
cox_table(
  data,
  time,
  event,
  factors = NULL,
  max_time = NULL,
  multivariable = c("all", "significant", "none"),
  significance = 0.05,
  ci_level = 0.95
)

logistic_table(
  data,
  outcome,
  factors = NULL,
  positive = NULL,
  multivariable = c("all", "significant", "none"),
  significance = 0.05,
  ci_level = 0.95
)
```

## Arguments

- data:

  Data frame holding the variables.

- time, event:

  Column names of the follow-up time and event indicator used by
  `cox_table()`.

- factors:

  Character vector of variable names to analyse. Defaults to every
  column of `data` except the time and event columns.

- max_time:

  Optional follow-up limit for `cox_table()`. Times beyond the limit are
  censored at the limit.

- multivariable:

  Which variables enter the multivariable model: `"all"` (default),
  `"significant"`, or `"none"`.

- significance:

  Significance level used by `multivariable = "significant"`.

- ci_level:

  Confidence level of the reported intervals.

- outcome:

  Column name of the binary outcome used by `logistic_table()`.

- positive:

  Value of `outcome` that represents the event.

## Value

A tibble with one row per term and the columns `term`, `variable`,
`level`, `n`, `n_event`, the univariable estimate (`hr_univariable` or
`or_univariable`) with interval and p-value, and the corresponding
multivariable columns.

## Details

The multivariable model includes all requested variables by default.
With `multivariable = "significant"` only variables whose univariable
likelihood-ratio test is below `significance` are included, which
reproduces the screening rule of the original helper functions.

## Examples

``` r
set.seed(42)
dat <- data.frame(
  time = rexp(60, 0.1),
  event = rbinom(60, 1, 0.7),
  age = rnorm(60, 60, 8),
  stage = factor(sample(c("II", "III"), 60, replace = TRUE))
)
cox_table(dat, time = "time", event = "event")
#> # A tibble: 2 × 13
#>   term   variable level     n n_event hr_univariable univ_ci_lower univ_ci_upper
#>   <chr>  <chr>    <chr> <int>   <int>          <dbl>         <dbl>         <dbl>
#> 1 age    age      NA       60      36          1.00          0.961          1.05
#> 2 stage… stage    III      60      36          0.968         0.495          1.89
#> # ℹ 5 more variables: univ_p <dbl>, hr_multivariable <dbl>,
#> #   multi_ci_lower <dbl>, multi_ci_upper <dbl>, multi_p <dbl>
set.seed(42)
dat <- data.frame(
  response = rbinom(60, 1, 0.4),
  age = rnorm(60, 60, 8),
  sex = factor(sample(c("F", "M"), 60, replace = TRUE))
)
logistic_table(dat, outcome = "response")
#> # A tibble: 2 × 13
#>   term  variable level     n n_event or_univariable univ_ci_lower univ_ci_upper
#>   <chr> <chr>    <chr> <int>   <int>          <dbl>         <dbl>         <dbl>
#> 1 age   age      NA       60      31          1.03          0.961          1.11
#> 2 sexM  sex      M        60      31          0.729         0.256          2.08
#> # ℹ 5 more variables: univ_p <dbl>, or_multivariable <dbl>,
#> #   multi_ci_lower <dbl>, multi_ci_upper <dbl>, multi_p <dbl>
```
