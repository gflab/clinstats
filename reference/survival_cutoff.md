# Survival cut point optimising the Youden index

Estimates Kaplan-Meier survival above and below each candidate cut point
of a marker at `predict_time` and returns the cut point that maximises
the Youden index (sensitivity + specificity - 1), the same criterion
that the long-unmaintained `survivalROC` package uses, computed here
directly from
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html).
Candidates are drawn from a quantile grid to keep the computation
bounded.

## Usage

``` r
survival_cutoff(time, event, marker, predict_time, grid_size = 200)
```

## Arguments

- time:

  Numeric follow-up times.

- event:

  Event indicator (0/1, logical, or a two-level value).

- marker:

  Numeric marker to split.

- predict_time:

  Time point at which sensitivity and specificity are evaluated.

- grid_size:

  Number of candidate cut points evaluated between the 1st and 99th
  percentile of the marker.

## Value

A tibble with one row: the optimal `cutoff`, the Youden index
(`youden`), the time-dependent true positive rate (`tp`) and false
positive rate (`fp`) at that cut point, and the number of observations.

## Examples

``` r
set.seed(1)
time <- rexp(100)
event <- rbinom(100, 1, 0.6)
marker <- rnorm(100)
survival_cutoff(time, event, marker, predict_time = 1)
#> # A tibble: 1 × 7
#>   cutoff youden    tp    fp predict_time     n n_event
#>    <dbl>  <dbl> <dbl> <dbl>        <dbl> <int>   <int>
#> 1   2.37  0.579     1 0.421            1   100      68
```
