# Oncotype DX colon cancer recurrence score

Computes the Oncotype DX colon cancer recurrence score from the
published gene panel. The formula and the panel, including the reference
genes, are unchanged from the original `calc_oncotypedx_crc()` helper.

## Usage

``` r
oncotype_crc(expression)
```

## Arguments

- expression:

  Numeric matrix with gene symbols in rows and samples in columns. Use
  [`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
  first if the matrix has samples in rows.

## Value

A tibble with one row per sample: the stromal, cell-cycle, and
individual gene scores, the reference score, the corrected scores, the
recurrence score, and the recurrence score category (`"Low"`,
`"Intermediate"`, or `"High"`).

## Examples

``` r
panel <- c(
  "BGN", "FAP", "INHBA", "MKI67", "MYC", "MYBL2", "GADD45B",
  "ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"
)
expr <- matrix(rnorm(12 * 3, 10), nrow = 12, dimnames = list(panel, paste0("s", 1:3)))
oncotype_crc(expr)
#> # A tibble: 3 × 11
#>   sample stroma cell_cycle individual reference corrected_stroma
#>   <chr>   <dbl>      <dbl>      <dbl>     <dbl>            <dbl>
#> 1 s1       10.2       9.80      11.1      10.2              9.96
#> 2 s2       11.1       9.57       8.78      9.64            11.5 
#> 3 s3       10.6      11.0        9.36      9.94            10.7 
#> # ℹ 5 more variables: corrected_cell_cycle <dbl>, corrected_individual <dbl>,
#> #   rs_score <dbl>, oncotype_score <dbl>, oncotype_class <fct>
```
