# Extract clinical variables from a legacy TCGA clinical table

Converts the wide clinical export distributed with legacy TCGA data
releases (column names of the form
`patient.age_at_initial_pathologic_diagnosis`) into an analysis-ready
tibble with follow-up times in days. Three column layouts are supported:
the generic layout, the oesophageal carcinoma layout (`"esca"`), and the
stomach adenocarcinoma layout (`"stad"`).

## Usage

``` r
tcga_clinical(clin, cancer = c("default", "esca", "stad"))
```

## Arguments

- clin:

  Data frame with one row per patient in the legacy TCGA layout.

- cancer:

  Column layout: `"default"`, `"esca"` (the historical `ESAC = TRUE`
  variant), or `"stad"`.

## Value

A tibble with the patient barcode, age, sex, stage, node status, T
stage, disease-free survival, relapse and death status with their times,
and, where available, chemotherapy and response. The `"esca"` layout
additionally returns the histological subtype and tumour location.
