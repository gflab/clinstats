# Clean an expression matrix and collapse duplicated gene columns

Prepares a gene expression matrix for downstream analysis: drops unnamed
and empty columns, expands columns with `" /// "` multi-mappings, maps
gene identifiers when requested, and keeps the column with the largest
median absolute deviation when several columns map to the same gene.

## Usage

``` r
clean_expression(expression, id_from = NULL, id_to = "SYMBOL")
```

## Arguments

- expression:

  Numeric matrix or data frame with samples in rows and gene identifiers
  in columns.

- id_from:

  Optional identifier type of the column names, for example
  `"ENTREZID"`. When supplied, identifiers are translated to `id_to`
  with
  [org.Hs.eg.db::org.Hs.eg.db](https://rdrr.io/pkg/org.Hs.eg.db/man/org.Hs.egBASE.html),
  which must be installed.

- id_to:

  Target identifier type used when `id_from` is supplied.

## Value

A numeric matrix with samples in rows and unique gene identifiers in
columns.

## Examples

``` r
expr <- matrix(
  c(1, 2, 3, 4, 5, 6, NA, 8),
  nrow = 2,
  dimnames = list(NULL, c("A", "A", "B", ""))
)
clean_expression(expr)
#>      A B
#> [1,] 1 5
#> [2,] 2 6
```
