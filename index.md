# clinstats

Clinical Research Statistics Toolkit. A compact set of statistical
utilities that recur in clinical and translational research, together
with curated colorectal cancer clinical datasets assembled from public
repositories and a Quarto analysis report template.

The package is the maintained successor of the archived `gaofenglib` and
`curatedClinicalData` repositories. Statistical routines were rewritten
with input validation, tidy output, and tests; the datasets are carried
over without changing the stored values.

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

All functions return tibbles, validate their inputs, and report clear
errors.

## Quick start

``` r

library(clinstats)

# Discriminatory performance of a marker
roc_summary(clin_crc$age, clin_crc$rfs.event)

# Univariable and multivariable Cox table
cox_table(
  clin_crc,
  time = "rfs.delay",
  event = "rfs.event",
  factors = c("sex", "age", "tnm.stage", "cms")
)

# Logistic table for a binary endpoint
logistic_table(
  clin_crc,
  outcome = "rfs.event",
  factors = c("sex", "age", "tnm.stage")
)
```

## Data

`clin_crc`, `clin_crc_gf`, and `clin_crc_cell` hold clinical annotations
for 2,816 to 2,868 colorectal cancer samples from GEO series GSE14333,
GSE17538, GSE18088, GSE26906, GSE31595, GSE33113, GSE37892, GSE39084,
GSE39582 (the PETACC-3 cohort), and TCGA clinical exports. All sources
are public and the data contain no patient identifiers. See
`data-raw/README.md` for provenance and checksums. Biological and
clinical conclusions drawn from these datasets remain the responsibility
of the analyst.

## Analysis report template

[`copy_report_template()`](https://gflab.github.io/clinstats/reference/copy_report_template.md)
writes a Quarto document with defaults for HTML and PDF output into the
current project.

``` r

copy_report_template()
```

## Migrating from gaofenglib and curatedClinicalData

The old packages are retired. The mapping below covers the exported
helpers; `init_R()` and `update_all()` were removed because packages
should not reinstall the user environment.

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

Behavioural notes for migrated code:

- Regression tables default to including all variables in the
  multivariable model. Use `multivariable = "significant"` to reproduce
  the previous univariable-screening rule.
- [`survival_cutoff()`](https://gflab.github.io/clinstats/reference/survival_cutoff.md)
  uses a Kaplan-Meier estimate of the time-dependent Youden index
  computed with the `survival` package, so the long-unmaintained
  `survivalROC` package is no longer a dependency.
- [`resample_cox()`](https://gflab.github.io/clinstats/reference/resample_cox.md)
  reports per-marker selection frequencies instead of a raw p-value
  matrix and follows the `future` parallel plan.
- [`oncotype_crc()`](https://gflab.github.io/clinstats/reference/oncotype_crc.md)
  expects genes in rows; transpose matrices that have samples in rows.

## License and citation

Apache License 2.0. See `CITATION.cff` for citation metadata. Copyright
2026 FengGao Lab contributors.

## Maintenance

Primary maintainer: Feng Gao (gflab). Bug reports and feature requests
are welcome through GitHub Issues.
