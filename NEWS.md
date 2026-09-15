# clinstats 0.1.0

* Initial release. Consolidates and modernises the statistical utilities and
  curated colorectal cancer clinical datasets previously distributed as the
  separate `gaofenglib` and `curatedClinicalData` repositories.
* Function redesign:
  * `roc_summary()` replaces `calc_logit()`, with tidy one-row output.
  * `odds_ratio()` replaces `calc_or()`.
  * `cox_table()` replaces `factor_analysis_cox()`.
  * `logistic_table()` replaces `factor_analysis_logit()`.
  * `evaluate_model()` replaces `eval_logit()`, with train/test evaluation and
    a DeLong comparison between the two areas under the curve.
  * `survival_cutoff()` replaces `calc_cutoff_survivalroc()` and computes the
    Kaplan-Meier Youden cut point directly with `survival::survfit()` over a
    bounded quantile grid, so the long-unmaintained `survivalROC` package is
    no longer a dependency.
  * `resample_cox()` replaces `calc_resamp_cox()`, reports per-gene selection
    frequencies, and defaults to the sequential `future` plan.
  * `clean_expression()` replaces `clean_dat()`; gene identifier mapping is
    now an optional `org.Hs.eg.db` dependency.
  * `oncotype_crc()` replaces `calc_oncotypedx_crc()`.
  * `tcga_clinical()` replaces `extract_tcga_clinical()`.
  * `survival_response()` replaces `get_survival()`.
* Removed environment-mutating helpers `init_R()` and `update_all()`.
* All exported functions return tibbles and validate their inputs with clear
  error messages.
* `cox_table()` and `logistic_table()` include every requested variable in
  the multivariable model by default; `multivariable = "significant"`
  reproduces the previous univariable-screening rule of the original
  helpers.
* Datasets `clin_crc`, `clin_crc_gf`, and `clin_crc_cell` are carried over
  byte-identical from `curatedClinicalData`, with provenance documentation.
* Added a Quarto analysis report template, `testthat` coverage, GitHub Actions
  checks, and a `pkgdown` site.
