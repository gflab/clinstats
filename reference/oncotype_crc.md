# Oncotype DX colon cancer recurrence score

Computes the Oncotype DX colon cancer recurrence score from the
published gene panel. The formula and the panel, including the reference
genes, are unchanged from the original `calc_oncotypedx_crc()` helper.

## Usage

``` r
oncotype_crc(expression, gene_axis = c("columns", "rows"))
```

## Arguments

- expression:

  Numeric matrix of gene expression values.

- gene_axis:

  Axis that holds the gene symbols. `"columns"` (the default) expects
  samples in rows, which is the orientation returned by
  [`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
  and the convention inherited from `gaofenglib`; the two functions
  therefore compose directly. Use `"rows"` for matrices with genes in
  rows and samples in columns.

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
samples_by_genes <- matrix(
  rnorm(3 * 12, 10),
  nrow = 3,
  dimnames = list(paste0("s", 1:3), panel)
)
oncotype_crc(samples_by_genes)
#> # A tibble: 3 × 11
#>   sample stroma cell_cycle individual reference corrected_stroma
#>   <chr>   <dbl>      <dbl>      <dbl>     <dbl>            <dbl>
#> 1 s1       9.66      10.7        8.78      9.80             9.85
#> 2 s2      10.1        9.63       9.53     10.1              9.94
#> 3 s3      10.6       10.7        9.38     10.7              9.91
#> # ℹ 5 more variables: corrected_cell_cycle <dbl>, corrected_individual <dbl>,
#> #   rs_score <dbl>, oncotype_score <dbl>, oncotype_class <fct>

# Genes in rows also work when declared explicitly.
oncotype_crc(t(samples_by_genes), gene_axis = "rows")
#> # A tibble: 3 × 11
#>   sample stroma cell_cycle individual reference corrected_stroma
#>   <chr>   <dbl>      <dbl>      <dbl>     <dbl>            <dbl>
#> 1 s1       9.66      10.7        8.78      9.80             9.85
#> 2 s2      10.1        9.63       9.53     10.1              9.94
#> 3 s3      10.6       10.7        9.38     10.7              9.91
#> # ℹ 5 more variables: corrected_cell_cycle <dbl>, corrected_individual <dbl>,
#> #   rs_score <dbl>, oncotype_score <dbl>, oncotype_class <fct>
```
