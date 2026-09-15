# Odds ratio for a marker with a binary outcome

Binarises a marker and reports the odds ratio for the event in the high
group relative to the low group. Continuous markers are split at the
Youden-optimal ROC threshold unless `cutoff` is supplied; two-level
markers and factors are used as they are.

## Usage

``` r
odds_ratio(x, outcome, positive = NULL, cutoff = NULL, ci_level = 0.95)
```

## Arguments

- x:

  Marker or factor to binarise. Continuous numeric values are split at
  the optimal threshold; numeric, logical, character, and factor inputs
  with two distinct values are used directly.

- outcome:

  Binary outcome. Logical vectors are coded `TRUE` as the event; numeric
  0/1 vectors use 1 as the event. Any other coding requires `positive`
  to name the event level.

- positive:

  Value of `outcome` that represents the event.

- cutoff:

  Optional threshold supplied by the user. Values greater than `cutoff`
  form the high group.

- ci_level:

  Confidence level for the AUC interval.

## Value

A tibble with one row: the threshold that was used, sample size, number
of events, group sizes, the odds ratio with a Wald confidence interval,
and the Wald p-value.

## Examples

``` r
set.seed(1)
marker <- c(rnorm(40, 0), rnorm(40, 1))
outcome <- rep(c(0, 1), each = 40)
odds_ratio(marker, outcome)
#> # A tibble: 1 × 9
#>   cutoff     n n_event n_high n_low odds_ratio ci_lower ci_upper   p_value
#>    <dbl> <int>   <int>  <int> <int>      <dbl>    <dbl>    <dbl>     <dbl>
#> 1  0.828    80      40     33    47       8.76     3.09     24.8 0.0000454
```
