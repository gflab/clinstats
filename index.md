# clinstats

Clinical research statistics, as tidy data. `clinstats` covers the
analyses that recur in clinical and translational papers —
discrimination, odds ratios, univariable and multivariable regression
tables, survival cut points, and expression-matrix preparation — and
returns a tibble from every function so results can be printed,
filtered, joined, and exported.

The package also ships curated colorectal cancer clinical datasets
assembled from public repositories, and a Quarto analysis report
template.

``` r

library(clinstats)

cox_table(
  clin_crc,
  time = "rfs.delay",
  event = "rfs.event",
  factors = c("sex", "age", "tnm.stage")
)
#> # A tibble: 3 × 13
#>   term      variable  level     n n_event hr_univariable univ_ci_lower
#>   <chr>     <chr>     <chr> <int>   <int>          <dbl>         <dbl>
#> 1 sexM      sex       M      1641     561          1.14          0.967
#> 2 age       age       <NA>   1640     561          0.999         0.993
#> 3 tnm.stage tnm.stage <NA>   1632     558          1.82          1.62
#> # ℹ 6 more variables: univ_ci_upper <dbl>, univ_p <dbl>,
#> #   hr_multivariable <dbl>, multi_ci_lower <dbl>, multi_ci_upper <dbl>,
#> #   multi_p <dbl>
```

## Why clinstats exists

Clinical statistics code tends to be written once per paper and copied
forward. The same ROC summary, the same two regression tables, the same
expression cleanup, reimplemented slightly differently each time — which
makes results hard to reproduce and easy to get subtly wrong.

`clinstats` collects those analyses into tested functions. Three
properties do the work:

- **A tibble from every function.** Results can be filtered, joined, and
  written out with the rest of a tidy pipeline, and printed in a report
  without reshaping. Regression tables put univariable and multivariable
  estimates side by side in one table.
- **Validated input.** Every function checks its arguments and reports
  what is wrong and how to fix it, instead of failing deep inside a
  model fit.
- **One expression-matrix orientation.** Samples in rows and genes in
  columns throughout, so the functions chain without manual transposes.

The statistical routines were rewritten from the retired `gaofenglib`
package, and the datasets carried over unchanged from
`curatedClinicalData`. See [Compatibility](#compatibility) for the
migration table.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("gflab/clinstats")
```

R 4.1 or later is required.

## What is inside

| Area | Functions |
|----|----|
| Discrimination | [`roc_summary()`](https://gflab.github.io/clinstats/reference/roc_summary.md), [`odds_ratio()`](https://gflab.github.io/clinstats/reference/odds_ratio.md) |
| Regression tables | [`cox_table()`](https://gflab.github.io/clinstats/reference/cox_table.md), [`logistic_table()`](https://gflab.github.io/clinstats/reference/cox_table.md), [`evaluate_model()`](https://gflab.github.io/clinstats/reference/evaluate_model.md) |
| Survival analysis | [`survival_cutoff()`](https://gflab.github.io/clinstats/reference/survival_cutoff.md), [`survival_response()`](https://gflab.github.io/clinstats/reference/survival_response.md), [`resample_cox()`](https://gflab.github.io/clinstats/reference/resample_cox.md) |
| Expression data | [`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md), [`oncotype_crc()`](https://gflab.github.io/clinstats/reference/oncotype_crc.md) |
| TCGA clinical data | [`tcga_clinical()`](https://gflab.github.io/clinstats/reference/tcga_clinical.md) |
| Reporting | [`copy_report_template()`](https://gflab.github.io/clinstats/reference/copy_report_template.md) |

## Quick start

### Discrimination and odds ratios

``` r

marker <- scale(clin_crc$age)

roc_summary(marker, clin_crc$rfs.event)
#> # A tibble: 1 × 12
#>       n n_event   auc auc_ci_lower auc_ci_upper direction threshold
#>   <int>   <int> <dbl>        <dbl>        <dbl> <chr>         <dbl>
#> 1  1698     574 0.514        0.485        0.542 >             0.143

odds_ratio(marker, clin_crc$rfs.event)
```

[`roc_summary()`](https://gflab.github.io/clinstats/reference/roc_summary.md)
returns the area under the curve, its confidence interval, the optimal
operating point, and the sensitivity and specificity at that point.
`direction = "auto"` lets `pROC` choose the orientation from the data,
which is what you want when the sign of a marker is not known in
advance.

### Regression tables

[`cox_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
and
[`logistic_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
fit every variable on its own, fit one multivariable model, and return
both sets of estimates in a single table. Numeric variables contribute
one row; factors contribute one row per level against the reference
level.

``` r

cox_table(
  clin_crc,
  time = "rfs.delay",
  event = "rfs.event",
  factors = c("sex", "age", "tnm.stage", "cms")
)
```

The multivariable model includes every requested variable by default.
Use `multivariable = "significant"` to keep only variables with a
univariable likelihood-ratio p-value below 0.05, or
`multivariable = "none"` for univariable estimates alone.

### Model evaluation

``` r

evaluate_model(train, test, outcome ~ marker1 + marker2)
```

Returns the training and test AUC with confidence intervals and a DeLong
comparison between them.

### Expression data

``` r

cleaned <- clean_expression(expression)     # samples in rows, genes in columns
oncotype_crc(cleaned)                       # composes directly
```

[`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
drops unusable columns, expands `" /// "` multi-mappings, optionally
maps identifiers with `org.Hs.eg.db`, and collapses duplicated gene
columns to the probe with the largest median absolute deviation.

## Data

`clin_crc`, `clin_crc_gf`, and `clin_crc_cell` hold clinical annotations
for 2,816 to 2,868 colorectal cancer samples across 25 columns:
demographics, TNM stage and T/N/M categories, tumour location and grade,
adjuvant chemotherapy, relapse-free, overall, and disease-free survival
with their times, MMR/CIMP/CIN status, TP53/KRAS/BRAF mutation status,
source dataset, and consensus molecular subtype.

Sources are the GEO series GSE14333, GSE17538, GSE18088, GSE26906,
GSE31595, GSE33113, GSE37892, GSE39084, and GSE39582 (the PETACC-3
cohort), plus TCGA clinical exports. Every source is public and the data
contain no patient identifiers.

``` r

dim(clin_crc)          # 2816 x 25
table(clin_crc$dataset)
```

[`?clin_crc`](https://gflab.github.io/clinstats/reference/clin_crc.md)
documents the column layout; `data-raw/README.md` records the source of
each cohort and the SHA-256 of each released object. Biological and
clinical conclusions drawn from these datasets remain the responsibility
of the analyst.

## Analysis report template

[`copy_report_template()`](https://gflab.github.io/clinstats/reference/copy_report_template.md)
writes a Quarto document into the current project, with settings for
HTML and PDF output and example chunks that use the bundled datasets.

``` r

copy_report_template()
```

## Compatibility

`clinstats` is the maintained successor of the archived `gaofenglib` and
`curatedClinicalData` repositories. Those packages still install, and
`clinstats` consolidates them so the lab maintains one package rather
than several overlapping ones.

| Old | New |
|----|----|
| `calc_logit()` | [`roc_summary()`](https://gflab.github.io/clinstats/reference/roc_summary.md) |
| `calc_or()` | [`odds_ratio()`](https://gflab.github.io/clinstats/reference/odds_ratio.md) |
| `factor_analysis_cox()` | [`cox_table()`](https://gflab.github.io/clinstats/reference/cox_table.md) |
| `factor_analysis_logit()` | [`logistic_table()`](https://gflab.github.io/clinstats/reference/cox_table.md) |
| `eval_logit()` | [`evaluate_model()`](https://gflab.github.io/clinstats/reference/evaluate_model.md) |
| `calc_cutoff_survivalroc()` | [`survival_cutoff()`](https://gflab.github.io/clinstats/reference/survival_cutoff.md) |
| `calc_resamp_cox()` | [`resample_cox()`](https://gflab.github.io/clinstats/reference/resample_cox.md) |
| `clean_dat()` | [`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md) |
| `calc_oncotypedx_crc()` | [`oncotype_crc()`](https://gflab.github.io/clinstats/reference/oncotype_crc.md) |
| `extract_tcga_clinical()` | [`tcga_clinical()`](https://gflab.github.io/clinstats/reference/tcga_clinical.md) |
| `get_survival()` | [`survival_response()`](https://gflab.github.io/clinstats/reference/survival_response.md) |
| `clin_crc`, `clin_crc_gf`, `clin_crc_cell` | same names |

Changes to expect when migrating:

- **Results are tibbles, not base data frames.** Positional indexing and
  [`rownames()`](https://rdrr.io/r/base/colnames.html) no longer apply.
- **Regression tables include every variable by default.** Use
  `multivariable = "significant"` to reproduce the previous
  univariable-screening rule.
- **[`survival_cutoff()`](https://gflab.github.io/clinstats/reference/survival_cutoff.md)
  computes the Kaplan-Meier Youden index** with the `survival` package,
  so the long-unmaintained `survivalROC` is no longer a dependency.
- **[`resample_cox()`](https://gflab.github.io/clinstats/reference/resample_cox.md)
  reports per-marker selection frequencies** instead of a raw p-value
  matrix, and follows the `future` parallel plan.
- **`init_R()` and `update_all()` were removed.** A library should not
  reinstall the user’s environment.

## Related packages

[gfplot](https://github.com/gflab/gfplot) produces the figures that go
with these analyses — Kaplan-Meier curves with hazard ratios, ROC
curves, risk score distributions — in the lab house style. The two
packages share the same datasets, so `clinstats` results feed directly
into `gfplot` figures.

## Citation

``` r

citation("clinstats")
```

See `CITATION.cff` for machine-readable metadata. Release `v0.1.2` is
available as a tag, so a specific state can be pinned.

## Getting help

Bug reports and feature requests are welcome through [GitHub
Issues](https://github.com/gflab/clinstats/issues). Reference
documentation for every function is on the [package
site](https://gflab.github.io/clinstats/).

## License

Apache License 2.0. See
[LICENSE.md](https://gflab.github.io/clinstats/LICENSE.md). Copyright
2026 FengGao Lab contributors.
