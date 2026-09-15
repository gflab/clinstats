#' Curated colorectal cancer clinical data
#'
#' Clinical information for colorectal cancer samples assembled from public
#' repositories: the GEO datasets GSE14333, GSE17538, GSE18088, GSE26906,
#' GSE31595, GSE33113, GSE37892, GSE39084, and GSE39582 (the PETACC-3
#' cohort), and TCGA colon and rectal cancer clinical exports. The data
#' contain research-level variables only and no patient identifiers. These
#' objects are carried over from the archived `curatedClinicalData`
#' repository without changing the stored values.
#'
#' @format A data frame with 25 columns:
#' \describe{
#'   \item{sample}{Sample or patient identifier.}
#'   \item{sex}{Sex (`"M"`/`"F"`).}
#'   \item{age}{Age at diagnosis in years.}
#'   \item{tnm.stage}{TNM stage as an ordered number.}
#'   \item{tnm.t, tnm.n, tnm.m}{TNM T, N, and M categories.}
#'   \item{lymphnodes}{Lymph node involvement.}
#'   \item{tumor.location}{Tumour location (`"L"`/`"R"`).}
#'   \item{grade}{Tumour grade.}
#'   \item{chemotherapy.adjuvant}{Adjuvant chemotherapy indicator.}
#'   \item{rfs.event, rfs.delay}{Relapse-free survival event and time
#'     (years).}
#'   \item{os.event, os.delay}{Overall survival event and time (years).}
#'   \item{dfs.event, dfs.delay}{Disease-free survival event and time
#'     (years).}
#'   \item{mmr.status, cimp.status, cin.status}{Molecular subtype
#'     indicators.}
#'   \item{tp53.mutation, kras.mutation, braf.mutation}{Mutation status.}
#'   \item{dataset}{Source dataset.}
#'   \item{cms}{Consensus molecular subtype label.}
#' }
#' @source GEO series GSE14333, GSE17538, GSE18088, GSE26906, GSE31595,
#'   GSE33113, GSE37892, GSE39084, and GSE39582; TCGA clinical exports.
#' @examples
#' dim(clin_crc)
#' table(clin_crc$dataset)
"clin_crc"

#' @rdname clin_crc
#' @format A data frame with 2,868 rows and 25 columns, in the same layout as
#'   [clin_crc]. This version adds the TCGA update of 2018.
"clin_crc_gf"

#' @rdname clin_crc
#' @format A data frame with 2,865 rows and 25 columns, in the same layout as
#'   [clin_crc]. This version carries the 2019 TCGA cell-line update.
"clin_crc_cell"
