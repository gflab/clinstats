# Getting started with clinstats

``` r

library(clinstats)
```

`clinstats` collects the statistical routines that recur in clinical
research workflows: discrimination summaries, odds ratios, regression
tables, survival cut points, and expression-matrix preparation. Every
function returns a tibble, so results can be printed for a report or
filtered, joined, and exported programmatically.

## Bundled data

The package ships a curated colorectal cancer clinical dataset assembled
from public GEO series and TCGA exports.

``` r

dim(clin_crc)
#> [1] 2816   25
table(clin_crc$dataset)
#> 
#> GSE14333 GSE17538 GSE18088 GSE26906 GSE31595 GSE33113 GSE37892 GSE39084 
#>      290      232       53       90       37       90      130       68 
#> GSE39582  PETACC3     TCGA 
#>      566      688      572
```

The three stored versions (`clin_crc`, `clin_crc_gf`, `clin_crc_cell`)
correspond to the 2018 and 2019 harmonisations of the original
repository. See the dataset help page for the column layout.

## Regression tables

[`logistic_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
and
[`cox_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
fit every variable on its own, fit one multivariable model, and return a
single table with both sets of estimates.

``` r

logistic_table(
  clin_crc,
  outcome = "rfs.event",
  factors = c("sex", "age", "tnm.stage", "cms")
)
#> # A tibble: 6 × 13
#>   term   variable level     n n_event or_univariable univ_ci_lower univ_ci_upper
#>   <chr>  <chr>    <chr> <int>   <int>          <dbl>         <dbl>         <dbl>
#> 1 sexM   sex      M      1699     574          1.11          0.911          1.36
#> 2 age    age      NA     1698     574          0.997         0.990          1.00
#> 3 tnm.s… tnm.sta… NA     1637     558          1.60          1.40           1.82
#> 4 cmsCM… cms      CMS2   1546     526          0.927         0.692          1.24
#> 5 cmsCM… cms      CMS3   1546     526          0.961         0.668          1.38
#> 6 cmsCM… cms      CMS4   1546     526          1.32          0.966          1.81
#> # ℹ 5 more variables: univ_p <dbl>, or_multivariable <dbl>,
#> #   multi_ci_lower <dbl>, multi_ci_upper <dbl>, multi_p <dbl>
```

``` r

cox_table(
  clin_crc,
  time = "rfs.delay",
  event = "rfs.event",
  factors = c("sex", "age", "tnm.stage", "cms")
)
#> # A tibble: 6 × 13
#>   term   variable level     n n_event hr_univariable univ_ci_lower univ_ci_upper
#>   <chr>  <chr>    <chr> <int>   <int>          <dbl>         <dbl>         <dbl>
#> 1 sexM   sex      M      1641     561          1.14          0.967          1.35
#> 2 age    age      NA     1640     561          0.999         0.993          1.01
#> 3 tnm.s… tnm.sta… NA     1632     558          1.82          1.62           2.04
#> 4 cmsCM… cms      CMS2   1493     514          0.868         0.681          1.11
#> 5 cmsCM… cms      CMS3   1493     514          0.877         0.650          1.18
#> 6 cmsCM… cms      CMS4   1493     514          1.30          1.01           1.68
#> # ℹ 5 more variables: univ_p <dbl>, hr_multivariable <dbl>,
#> #   multi_ci_lower <dbl>, multi_ci_upper <dbl>, multi_p <dbl>
```

The multivariable model includes all requested variables by default.
Pass `multivariable = "significant"` to keep the historical screening
rule in which only variables with a univariable likelihood-ratio p-value
below 0.05 enter the model.

## Marker evaluation

[`roc_summary()`](https://gflab.github.io/clinstats/reference/roc_summary.md)
reports the area under the curve together with the optimal operating
point, and
[`odds_ratio()`](https://gflab.github.io/clinstats/reference/odds_ratio.md)
reports the odds ratio for the high group at the Youden-optimal
threshold.

``` r

marker <- scale(clin_crc$age)
roc_summary(marker, clin_crc$rfs.event)
#> # A tibble: 1 × 12
#>       n n_event   auc auc_ci_lower auc_ci_upper direction threshold sensitivity
#>   <int>   <int> <dbl>        <dbl>        <dbl> <chr>         <dbl>       <dbl>
#> 1  1698     574 0.514        0.485        0.542 >             0.143       0.493
#> # ℹ 4 more variables: specificity <dbl>, ppv <dbl>, npv <dbl>, accuracy <dbl>
odds_ratio(marker, clin_crc$rfs.event)
#> # A tibble: 1 × 9
#>   cutoff     n n_event n_high n_low odds_ratio ci_lower ci_upper p_value
#>    <dbl> <int>   <int>  <int> <int>      <dbl>    <dbl>    <dbl>   <dbl>
#> 1  0.143  1698     574    913   785      0.830    0.678     1.02  0.0697
```

For a survival endpoint,
[`survival_cutoff()`](https://gflab.github.io/clinstats/reference/survival_cutoff.md)
estimates the Kaplan-Meier Youden index over a grid of candidate cut
points and returns the best one.

``` r

survival_cutoff(
  time = clin_crc$rfs.delay,
  event = clin_crc$rfs.event,
  marker = marker,
  predict_time = 3
)
#> # A tibble: 1 × 7
#>   cutoff  youden     tp     fp predict_time     n n_event
#>    <dbl>   <dbl>  <dbl>  <dbl>        <dbl> <int>   <int>
#> 1   1.49 -0.0181 0.0284 0.0465            3  1640     561
```

## Expression data

[`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
prepares an expression matrix: it drops unusable columns, expands
`" /// "` multi-mappings, optionally maps identifiers with
`org.Hs.eg.db`, and collapses duplicated gene columns to the probe with
the largest median absolute deviation.

``` r

expr <- matrix(
  rnorm(6 * 4),
  nrow = 6,
  dimnames = list(NULL, c("A", "A", "B", "C /// D"))
)
clean_expression(expr)
#>                 A           B          C          D
#> [1,] -1.400043517  2.06502490  0.5429963  0.5429963
#> [2,]  0.255317055 -1.63098940 -0.9140748 -0.9140748
#> [3,] -2.437263611  0.51242695  0.4681544  0.4681544
#> [4,] -0.005571287 -1.86301149  0.3629513  0.3629513
#> [5,]  0.621552721 -0.52201251 -1.3045435 -1.3045435
#> [6,]  1.148411606 -0.05260191  0.7377763  0.7377763
```

[`oncotype_crc()`](https://gflab.github.io/clinstats/reference/oncotype_crc.md)
computes the Oncotype DX colon cancer recurrence score from the
published panel. It expects the same orientation as
[`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
— samples in rows and genes in columns — so the two chain directly:

``` r

panel <- c(
  "BGN", "FAP", "INHBA", "MKI67", "MYC", "MYBL2", "GADD45B",
  "ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"
)
set.seed(1)
expression <- matrix(
  rnorm(5 * length(panel), 10),
  nrow = 5,
  dimnames = list(paste0("s", 1:5), panel)
)
oncotype_crc(clean_expression(expression))
#> # A tibble: 5 × 11
#>   sample stroma cell_cycle individual reference corrected_stroma
#>   <chr>   <dbl>      <dbl>      <dbl>     <dbl>            <dbl>
#> 1 s1      10.0       10.3       11.4      10.2              9.80
#> 2 s2      10.4       10.2        9.90      9.75            10.6 
#> 3 s3       9.76       9.85      10.4      10.1              9.62
#> 4 s4       9.99       9.45       9.95     10.2              9.79
#> 5 s5      10.4       10.5        8.62     10.5              9.93
#> # ℹ 5 more variables: corrected_cell_cycle <dbl>, corrected_individual <dbl>,
#> #   rs_score <dbl>, oncotype_score <dbl>, oncotype_class <fct>
```

Matrices with genes in rows are accepted by declaring the axis:
`oncotype_crc(t(expression), gene_axis = "rows")`.

## Reproducibility

The bundled datasets are byte-identical to the archived
`curatedClinicalData` objects, with provenance and checksums recorded in
`data-raw/README.md`. Statistical conclusions drawn from them remain the
responsibility of the analyst.
