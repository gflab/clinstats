# Changelog

## clinstats 0.1.1

- One orientation for every expression-matrix function: **samples in
  rows, genes in columns**, which is the convention inherited from
  `gaofenglib`. In 0.1.0
  [`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
  returned samples in rows while
  [`oncotype_crc()`](https://gflab.github.io/clinstats/reference/oncotype_crc.md)
  required genes in rows, so the two could not be chained without a
  manual transpose.
- [`oncotype_crc()`](https://gflab.github.io/clinstats/reference/oncotype_crc.md)
  gains a `gene_axis` argument. The default `"columns"` matches
  [`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
  and composes with it directly; `"rows"` keeps the previous orientation
  available. Declaring the wrong axis now reports which argument to use
  instead of listing twelve missing genes.

## clinstats 0.1.0

- Initial release. Consolidates and modernises the statistical utilities
  and curated colorectal cancer clinical datasets previously distributed
  as the separate `gaofenglib` and `curatedClinicalData` repositories.
- Function redesign:
  - [`roc_summary()`](https://gflab.github.io/clinstats/reference/roc_summary.md)
    replaces `calc_logit()`, with tidy one-row output.
  - [`odds_ratio()`](https://gflab.github.io/clinstats/reference/odds_ratio.md)
    replaces `calc_or()`.
  - [`cox_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
    replaces `factor_analysis_cox()`.
  - [`logistic_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
    replaces `factor_analysis_logit()`.
  - [`evaluate_model()`](https://gflab.github.io/clinstats/reference/evaluate_model.md)
    replaces `eval_logit()`, with train/test evaluation and a DeLong
    comparison between the two areas under the curve.
  - [`survival_cutoff()`](https://gflab.github.io/clinstats/reference/survival_cutoff.md)
    replaces `calc_cutoff_survivalroc()` and computes the Kaplan-Meier
    Youden cut point directly with
    [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
    over a bounded quantile grid, so the long-unmaintained `survivalROC`
    package is no longer a dependency.
  - [`resample_cox()`](https://gflab.github.io/clinstats/reference/resample_cox.md)
    replaces `calc_resamp_cox()`, reports per-gene selection
    frequencies, and defaults to the sequential `future` plan.
  - [`clean_expression()`](https://gflab.github.io/clinstats/reference/clean_expression.md)
    replaces `clean_dat()`; gene identifier mapping is now an optional
    `org.Hs.eg.db` dependency.
  - [`oncotype_crc()`](https://gflab.github.io/clinstats/reference/oncotype_crc.md)
    replaces `calc_oncotypedx_crc()`.
  - [`tcga_clinical()`](https://gflab.github.io/clinstats/reference/tcga_clinical.md)
    replaces `extract_tcga_clinical()`.
  - [`survival_response()`](https://gflab.github.io/clinstats/reference/survival_response.md)
    replaces `get_survival()`.
- Removed environment-mutating helpers `init_R()` and `update_all()`.
- All exported functions return tibbles and validate their inputs with
  clear error messages.
- [`cox_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
  and
  [`logistic_table()`](https://gflab.github.io/clinstats/reference/cox_table.md)
  include every requested variable in the multivariable model by
  default; `multivariable = "significant"` reproduces the previous
  univariable-screening rule of the original helpers.
- Datasets `clin_crc`, `clin_crc_gf`, and `clin_crc_cell` are carried
  over byte-identical from `curatedClinicalData`, with provenance
  documentation.
- Added a Quarto analysis report template, `testthat` coverage, GitHub
  Actions checks, and a `pkgdown` site.
