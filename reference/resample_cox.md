# Resampling-based Cox screening of candidate markers

Repeatedly draws subsamples of the samples, fits a univariable Cox model
for every candidate column, and reports how often each candidate reached
`p_threshold`. This is the resampling screening previously provided by
`calc_resamp_cox()`, reported as per-marker selection frequencies
instead of a raw p-value matrix.

## Usage

``` r
resample_cox(
  expression,
  time,
  event,
  ratio = 0.8,
  times = 1000,
  p_threshold = 0.05,
  seed = NULL
)
```

## Arguments

- expression:

  Numeric matrix or data frame with samples in rows and candidate
  markers in columns.

- time, event:

  Follow-up time and event indicator, one value per sample.

- ratio:

  Fraction of samples drawn in each resample.

- times:

  Number of resamples.

- p_threshold:

  Significance threshold that counts as a selection.

- seed:

  Optional integer seed for reproducible resampling.

## Value

A tibble with one row per candidate marker: the number of valid
resamples, the number of resamples in which the marker was significant,
the selection fraction, and the median p-value.

## Details

The function uses the `future` framework. Set a parallel plan before
calling it to use several workers, for example
`future::plan(future::multisession)`.

## Examples

``` r
set.seed(1)
expr <- matrix(rnorm(60 * 5), nrow = 60, dimnames = list(NULL, paste0("g", 1:5)))
time <- rexp(60)
event <- rbinom(60, 1, 0.7)
resample_cox(expr, time, event, times = 20, seed = 1)
#> # A tibble: 5 × 5
#>   gene  n_valid n_significant fraction median_p
#>   <chr>   <dbl>         <dbl>    <dbl>    <dbl>
#> 1 g1         20             1     0.05    0.140
#> 2 g2         20             0     0       0.661
#> 3 g3         20             0     0       0.794
#> 4 g4         20             0     0       0.297
#> 5 g5         20             0     0       0.715
```
