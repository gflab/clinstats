# Data provenance

`clin_crc`, `clin_crc_gf`, and `clin_crc_cell` were assembled in the
`curatedClinicalData` repository between 2018 and 2019 and are carried over
here without changing the stored values. This directory documents where the
data come from; the public GEO series can be re-downloaded with `GEOquery`,
and the TCGA clinical exports with `TCGAbiolinks`.

## Sources

* GEO series: GSE14333, GSE17538, GSE18088, GSE26906, GSE31595, GSE33113,
  GSE37892, GSE39084, GSE39582 (the PETACC-3 cohort).
* TCGA colon and rectal cancer clinical exports (2018 version in
  `clin_crc_gf`, 2019 cell-line update in `clin_crc_cell`).
* The `data-raw/crc.cit.clin.RData` intermediate in the original repository
  corresponds to the CIT (Cartes d'Identite des Tumeurs) cohort distributed
  with GSE39582.

All variables are research-level clinical annotations. No patient names,
record numbers, or other direct identifiers are present, and the lab policy
of not uploading patient-level data to GitHub is respected: every source is
a public repository.

## Checksums of the released objects

```
1920ae4e0f38fb879997de8d3d670ba81e6a75c27ae741b45b02310ed735dfe5  clin_crc_cell.rda
99b79d47fcde9516ba0c6f0b74367b3cd0183015615f6601b57e669b05386712  clin_crc_gf.rda
fdc436fdd03192a9eda954d18e961d5feca29e01761301bd8f5022d5ff91841d  clin_crc.rda
```

## Rebuilding

The original assembly script is preserved for reference; it documents the
mapping from GEO and TCGA exports to the harmonised columns and is not
required to use the package. Reproducing the exact objects requires the
2018/2019 versions of the source exports.
