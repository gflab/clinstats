test_that("clean_expression drops empty columns and keeps the best duplicate", {
  expr <- matrix(
    c(1, 2, 3, 0, 5, 10, NA, 8, 9, 10, 11, 12),
    nrow = 3,
    dimnames = list(NULL, c("A", "A", "B", ""))
  )
  result <- clean_expression(expr)
  expect_equal(colnames(result), c("A", "B"))
  expect_equal(result[, "A"], expr[, 2])
  expect_equal(result[, "B"], expr[, 3])
})

test_that("clean_expression expands multi-mapping columns", {
  expr <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(NULL, c("100 /// 200", "300"))
  )
  result <- clean_expression(expr)
  expect_equal(colnames(result), c("100", "200", "300"))
  expect_equal(result[, "100"], result[, "200"])
})

test_that("oncotype_crc follows the published formula", {
  panel <- c(
    "BGN", "FAP", "INHBA", "MKI67", "MYC", "MYBL2", "GADD45B",
    "ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"
  )
  set.seed(6)
  samples <- paste0("s", 1:5)
  expr <- matrix(
    rnorm(5 * 12, 10),
    nrow = 5,
    dimnames = list(samples, panel)
  )
  result <- oncotype_crc(expr)

  expect_equal(result$sample, samples)

  # Reference values are computed on the transposed matrix so the test does
  # not reuse the orientation the function normalises internally.
  genes_by_samples <- t(expr)
  manual_reference <- colMeans(
    genes_by_samples[c("ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"), ]
  )
  manual_stroma <- colMeans(genes_by_samples[c("BGN", "FAP", "INHBA"), ])
  manual_cell_cycle <- colMeans(genes_by_samples[c("MKI67", "MYC", "MYBL2"), ])
  manual_rs <- 0.15 * (manual_stroma - manual_reference + 10) -
    0.3 * (manual_cell_cycle - manual_reference + 10) +
    0.15 * (genes_by_samples["GADD45B", ] - manual_reference + 10)
  expect_equal(
    result$oncotype_score,
    unname(44 * (manual_rs + 0.82)),
    ignore_attr = TRUE
  )
  expect_true(all(levels(result$oncotype_class) == c("Low", "Intermediate", "High")))
})

test_that("clean_expression output feeds oncotype_crc without transposing", {
  panel <- c(
    "BGN", "FAP", "INHBA", "MKI67", "MYC", "MYBL2", "GADD45B",
    "ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"
  )
  set.seed(7)
  samples <- paste0("s", 1:6)
  expr <- matrix(
    rnorm(6 * length(panel), 10),
    nrow = 6,
    dimnames = list(samples, panel)
  )
  cleaned <- clean_expression(expr)
  result <- oncotype_crc(cleaned)

  expect_equal(nrow(result), 6L)
  expect_equal(result$sample, samples)
  expect_equal(
    result$oncotype_score,
    oncotype_crc(expr)$oncotype_score,
    ignore_attr = TRUE
  )
})

test_that("oncotype_crc accepts genes in rows and reports a useful hint", {
  panel <- c(
    "BGN", "FAP", "INHBA", "MKI67", "MYC", "MYBL2", "GADD45B",
    "ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"
  )
  set.seed(8)
  genes_by_samples <- matrix(
    rnorm(length(panel) * 4, 10),
    nrow = length(panel),
    dimnames = list(panel, paste0("s", 1:4))
  )
  samples_by_genes <- t(genes_by_samples)

  expect_equal(
    oncotype_crc(genes_by_samples, gene_axis = "rows")$oncotype_score,
    oncotype_crc(samples_by_genes)$oncotype_score,
    ignore_attr = TRUE
  )

  # Declaring the wrong axis points at the correct argument instead of
  # reporting twelve missing genes.
  expect_error(
    oncotype_crc(genes_by_samples),
    'gene_axis = "rows"'
  )
  expect_error(
    oncotype_crc(samples_by_genes, gene_axis = "rows"),
    'gene_axis = "columns"'
  )
})

test_that("oncotype_crc reports genuinely missing panel genes", {
  panel <- c(
    "BGN", "FAP", "INHBA", "MKI67", "MYC", "MYBL2", "GADD45B",
    "ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"
  )
  set.seed(9)
  expr <- matrix(
    rnorm(3 * 12, 10),
    nrow = 3,
    dimnames = list(paste0("s", 1:3), panel)
  )
  trimmed <- expr[, -1, drop = FALSE]

  expect_error(oncotype_crc(trimmed), "missing panel genes")
})

test_that("expression matrices keep samples in rows throughout", {
  # One orientation for the package: clean_expression() and oncotype_crc()
  # both treat rows as samples and columns as genes, and neither transposes
  # the sample axis behind the caller's back.
  expr <- matrix(
    rnorm(4 * 3),
    nrow = 4,
    dimnames = list(paste0("s", 1:4), c("A", "B", "C"))
  )
  result <- clean_expression(expr)
  expect_equal(rownames(result), rownames(expr))
  expect_equal(colnames(result), c("A", "B", "C"))
})
