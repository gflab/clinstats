#' Clean an expression matrix and collapse duplicated gene columns
#'
#' Prepares a gene expression matrix for downstream analysis: drops unnamed
#' and empty columns, expands columns with `" /// "` multi-mappings, maps
#' gene identifiers when requested, and keeps the column with the largest
#' median absolute deviation when several columns map to the same gene.
#'
#' @param expression Numeric matrix or data frame with samples in rows and
#'   gene identifiers in columns.
#' @param id_from Optional identifier type of the column names, for example
#'   `"ENTREZID"`. When supplied, identifiers are translated to `id_to` with
#'   [org.Hs.eg.db::org.Hs.eg.db], which must be installed.
#' @param id_to Target identifier type used when `id_from` is supplied.
#'
#' @return A numeric matrix with samples in rows and unique gene identifiers
#'   in columns.
#'
#' @export
#' @examples
#' expr <- matrix(
#'   c(1, 2, 3, 4, 5, 6, NA, 8),
#'   nrow = 2,
#'   dimnames = list(NULL, c("A", "A", "B", ""))
#' )
#' clean_expression(expr)
clean_expression <- function(expression, id_from = NULL, id_to = "SYMBOL") {
  x <- as.matrix(expression)
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg expression} must be numeric.")
  }
  keep <- !is.na(colnames(x)) & nzchar(colnames(x))
  x <- x[, keep, drop = FALSE]
  keep <- colSums(!is.na(x)) > 0
  x <- x[, keep, drop = FALSE]
  if (ncol(x) == 0) {
    cli::cli_abort("No usable columns remain in {.arg expression}.")
  }

  if (any(grepl("///", colnames(x), fixed = TRUE))) {
    ids <- strsplit(colnames(x), " /// ", fixed = TRUE)
    x <- do.call(cbind, lapply(seq_along(ids), function(i) {
      out <- matrix(
        rep(x[, i], times = length(ids[[i]])),
        ncol = length(ids[[i]])
      )
      colnames(out) <- trimws(ids[[i]])
      out
    }))
    keep <- colSums(!is.na(x)) > 0
    x <- x[, keep, drop = FALSE]
  }

  if (!is.null(id_from) && id_from != id_to) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      cli::cli_abort(c(
        "Package {.pkg org.Hs.eg.db} is required to map identifiers.",
        i = "Install it with {.code BiocManager::install(\"org.Hs.eg.db\")}."
      ))
    }
    mapped <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = colnames(x),
      keytype = id_from,
      column = id_to,
      multiVals = "first"
    ))
    colnames(x) <- as.character(mapped)
    keep <- !is.na(colnames(x)) & nzchar(colnames(x))
    x <- x[, keep, drop = FALSE]
  }

  mad_col <- suppressWarnings(apply(x, 2, stats::mad, na.rm = TRUE))
  mad_col[is.na(mad_col)] <- -Inf
  order_index <- order(colnames(x), -mad_col)
  x <- x[, order_index, drop = FALSE]
  x[, !duplicated(colnames(x)), drop = FALSE]
}

#' Oncotype DX colon cancer recurrence score
#'
#' Computes the Oncotype DX colon cancer recurrence score from the published
#' gene panel. The formula and the panel, including the reference genes, are
#' unchanged from the original `calc_oncotypedx_crc()` helper.
#'
#' @param expression Numeric matrix of gene expression values.
#' @param gene_axis Axis that holds the gene symbols. `"columns"` (the
#'   default) expects samples in rows, which is the orientation returned by
#'   [clean_expression()] and the convention inherited from `gaofenglib`;
#'   the two functions therefore compose directly. Use `"rows"` for matrices
#'   with genes in rows and samples in columns.
#'
#' @return A tibble with one row per sample: the stromal, cell-cycle, and
#'   individual gene scores, the reference score, the corrected scores, the
#'   recurrence score, and the recurrence score category (`"Low"`,
#'   `"Intermediate"`, or `"High"`).
#'
#' @export
#' @examples
#' panel <- c(
#'   "BGN", "FAP", "INHBA", "MKI67", "MYC", "MYBL2", "GADD45B",
#'   "ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"
#' )
#' samples_by_genes <- matrix(
#'   rnorm(3 * 12, 10),
#'   nrow = 3,
#'   dimnames = list(paste0("s", 1:3), panel)
#' )
#' oncotype_crc(samples_by_genes)
#'
#' # Genes in rows also work when declared explicitly.
#' oncotype_crc(t(samples_by_genes), gene_axis = "rows")
oncotype_crc <- function(expression, gene_axis = c("columns", "rows")) {
  gene_axis <- match.arg(gene_axis)
  x <- as.matrix(expression)
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg expression} must be numeric.")
  }

  stroma <- c("BGN", "FAP", "INHBA")
  cell_cycle <- c("MKI67", "MYC", "MYBL2")
  individual <- "GADD45B"
  reference <- c("ATP5E", "GPX1", "PGK1", "VDAC2", "UBB")
  panel <- c(stroma, cell_cycle, individual, reference)

  probes <- if (gene_axis == "columns") colnames(x) else rownames(x)
  missing <- setdiff(panel, probes)
  if (length(missing) > 0) {
    other <- if (gene_axis == "columns") rownames(x) else colnames(x)
    if (!is.null(other) && all(panel %in% other)) {
      alternative <- if (gene_axis == "columns") "rows" else "columns"
      cli::cli_abort(c(
        "The panel genes are on the {alternative} axis of {.arg expression}.",
        i = 'Pass {.code gene_axis = "{alternative}"}, or transpose the matrix.'
      ))
    }
    cli::cli_abort(c(
      "{.arg expression} is missing panel genes {.val {missing}}.",
      i = "Panel genes must be on the axis named by {.arg gene_axis}."
    ))
  }

  # Normalise once to genes in rows so the score is computed in one
  # orientation regardless of how the caller supplied the matrix.
  if (gene_axis == "columns") {
    sample_ids <- rownames(x)
    x <- t(x)
  } else {
    sample_ids <- colnames(x)
  }
  if (is.null(sample_ids)) {
    sample_ids <- as.character(seq_len(ncol(x)))
  }

  mean_score <- function(genes) {
    colMeans(x[genes, , drop = FALSE], na.rm = TRUE)
  }
  score_stroma <- mean_score(stroma)
  score_cell_cycle <- mean_score(cell_cycle)
  score_individual <- x[individual, ]
  score_reference <- mean_score(reference)
  corrected_stroma <- score_stroma - score_reference + 10
  corrected_cell_cycle <- score_cell_cycle - score_reference + 10
  corrected_individual <- score_individual - score_reference + 10
  rs_score <- 0.15 * corrected_stroma -
    0.3 * corrected_cell_cycle +
    0.15 * corrected_individual
  oncotype_score <- 44 * (rs_score + 0.82)

  tibble::tibble(
    sample = sample_ids,
    stroma = score_stroma,
    cell_cycle = score_cell_cycle,
    individual = score_individual,
    reference = score_reference,
    corrected_stroma = corrected_stroma,
    corrected_cell_cycle = corrected_cell_cycle,
    corrected_individual = corrected_individual,
    rs_score = rs_score,
    oncotype_score = oncotype_score,
    oncotype_class = cut(
      oncotype_score,
      breaks = c(-Inf, 30, 41, Inf),
      labels = c("Low", "Intermediate", "High")
    )
  )
}
