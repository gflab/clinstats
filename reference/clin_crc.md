# Curated colorectal cancer clinical data

Clinical information for colorectal cancer samples assembled from public
repositories: the GEO datasets GSE14333, GSE17538, GSE18088, GSE26906,
GSE31595, GSE33113, GSE37892, GSE39084, and GSE39582 (the PETACC-3
cohort), and TCGA colon and rectal cancer clinical exports. The data
contain research-level variables only and no patient identifiers. These
objects are carried over from the archived `curatedClinicalData`
repository without changing the stored values.

## Usage

``` r
clin_crc

clin_crc_gf

clin_crc_cell
```

## Format

A data frame with 25 columns:

- sample:

  Sample or patient identifier.

- sex:

  Sex (`"M"`/`"F"`).

- age:

  Age at diagnosis in years.

- tnm.stage:

  TNM stage as an ordered number.

- tnm.t, tnm.n, tnm.m:

  TNM T, N, and M categories.

- lymphnodes:

  Lymph node involvement.

- tumor.location:

  Tumour location (`"L"`/`"R"`).

- grade:

  Tumour grade.

- chemotherapy.adjuvant:

  Adjuvant chemotherapy indicator.

- rfs.event, rfs.delay:

  Relapse-free survival event and time (years).

- os.event, os.delay:

  Overall survival event and time (years).

- dfs.event, dfs.delay:

  Disease-free survival event and time (years).

- mmr.status, cimp.status, cin.status:

  Molecular subtype indicators.

- tp53.mutation, kras.mutation, braf.mutation:

  Mutation status.

- dataset:

  Source dataset.

- cms:

  Consensus molecular subtype label.

A data frame with 2,868 rows and 25 columns, in the same layout as
clin_crc. This version adds the TCGA update of 2018.

A data frame with 2,865 rows and 25 columns, in the same layout as
clin_crc. This version carries the 2019 TCGA cell-line update.

## Source

GEO series GSE14333, GSE17538, GSE18088, GSE26906, GSE31595, GSE33113,
GSE37892, GSE39084, and GSE39582; TCGA clinical exports.

## Examples

``` r
dim(clin_crc)
#> [1] 2816   25
table(clin_crc$dataset)
#> 
#> GSE14333 GSE17538 GSE18088 GSE26906 GSE31595 GSE33113 GSE37892 GSE39084 
#>      290      232       53       90       37       90      130       68 
#> GSE39582  PETACC3     TCGA 
#>      566      688      572 
```
