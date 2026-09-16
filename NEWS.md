# clinstats 0.1.2

* `evaluate_model()` reported `1 - AUC`. The ROC curves were built with
  `direction = ">"`, but pROC defines `"<"` as "observations are positive
  when they are greater than or equal to the threshold", which is the
  convention a predicted probability follows. A model that ranked the data
  well was therefore reported as worse than chance (a true AUC of 0.78 was
  reported as 0.22) and the DeLong comparison tested the inverted hypothesis.
  The direction is now `"<"`.
* Added a semantic test for this: a marker that is higher in the event group
  must produce an AUC above 0.5. The previous test compared the result
  against a `pROC::roc()` call that used the same inverted direction, so it
  could not detect the defect.

# clinstats 0.1.1

* One orientation for every expression-matrix function: **samples in rows,
  genes in columns**, which is the convention inherited from `gaofenglib`.
  In 0.1.0 `clean_expression()` returned samples in rows while
  `oncotype_crc()` required genes in rows, so the two could not be chained
  without a manual transpose.
* `oncotype_crc()` gains a `gene_axis` argument. The default `"columns"`
  matches `clean_expression()` and composes with it directly; `"rows"` keeps
  the previous orientation available. Declaring the wrong axis now reports
  which argument to use instead of listing twelve missing genes.

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
