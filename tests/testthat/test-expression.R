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
  expr <- matrix(
    rnorm(12 * 5, 10),
    nrow = 12,
    dimnames = list(panel, paste0("s", 1:5))
  )
  result <- oncotype_crc(expr)

  reference <- mean(expr["ATP5E", ] + expr["GPX1", ] + expr["PGK1", ] +
    expr["VDAC2", ] + expr["UBB", ]) / 5
  manual_reference <- colMeans(expr[c("ATP5E", "GPX1", "PGK1", "VDAC2", "UBB"), ])
  manual_stroma <- colMeans(expr[c("BGN", "FAP", "INHBA"), ])
  manual_cell_cycle <- colMeans(expr[c("MKI67", "MYC", "MYBL2"), ])
  manual_rs <- 0.15 * (manual_stroma - manual_reference + 10) -
    0.3 * (manual_cell_cycle - manual_reference + 10) +
    0.15 * (expr["GADD45B", ] - manual_reference + 10)
  expect_equal(
    result$oncotype_score,
    unname(44 * (manual_rs + 0.82)),
    ignore_attr = TRUE
  )
  expect_true(all(levels(result$oncotype_class) == c("Low", "Intermediate", "High")))

  expect_error(oncotype_crc(expr[-1, ]), "missing panel genes")
})
