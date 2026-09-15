# Summarise the diagnostic performance of a numeric marker

Computes the area under the receiver operating characteristic curve
(AUC) with a confidence interval and the optimal operating point of a
numeric marker against a binary outcome. A tibble with one row is
returned, so the result can be combined across markers or subgroups.

## Usage

``` r
roc_summary(
  marker,
  outcome,
  positive = NULL,
  direction = "auto",
  best_method = c("youden", "closest.topleft"),
  ci_level = 0.95,
  power = FALSE,
  alternative = c("one.sided", "two.sided")
)
```

## Arguments

- marker:

  Numeric vector of marker values.

- outcome:

  Binary outcome. Logical vectors are coded `TRUE` as the event; numeric
  0/1 vectors use 1 as the event. Any other coding requires `positive`
  to name the event level.

- positive:

  Value of `outcome` that represents the event.

- direction:

  Direction of the comparison passed on to
  [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html). `"auto"` lets
  pROC choose the direction from the data; `"<"` and `">"` fix it
  explicitly.

- best_method:

  Method used to choose the optimal operating point, `"youden"`
  (default) or `"closest.topleft"`.

- ci_level:

  Confidence level for the AUC interval.

- power:

  If `TRUE`, add the power of the AUC test as reported by
  [`pROC::power.roc.test()`](https://rdrr.io/pkg/pROC/man/power.roc.test.html).
  The value may be `NA` when power cannot be estimated for the sample.

- alternative:

  Alternative hypothesis used for the power estimate.

## Value

A tibble with one row: sample size, number of events, AUC and its
confidence interval, the optimal threshold, and the sensitivity,
specificity, positive and negative predictive values, and accuracy at
that threshold.

## See also

[`odds_ratio()`](https://gflab.github.io/clinstats/reference/odds_ratio.md)
for the odds ratio at the optimal threshold.

## Examples

``` r
set.seed(1)
marker <- c(rnorm(40, 0), rnorm(40, 1))
outcome <- rep(c(0, 1), each = 40)
roc_summary(marker, outcome)
#> # A tibble: 1 × 12
#>       n n_event   auc auc_ci_lower auc_ci_upper direction threshold sensitivity
#>   <int>   <int> <dbl>        <dbl>        <dbl> <chr>         <dbl>       <dbl>
#> 1    80      40 0.789        0.692        0.887 <             0.828        0.65
#> # ℹ 4 more variables: specificity <dbl>, ppv <dbl>, npv <dbl>, accuracy <dbl>
```
