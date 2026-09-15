test_that("curated datasets are intact", {
  expect_equal(dim(clin_crc), c(2816L, 25L))
  expect_equal(dim(clin_crc_gf), c(2868L, 25L))
  expect_equal(dim(clin_crc_cell), c(2865L, 25L))
  expect_equal(
    colnames(clin_crc),
    c(
      "sample", "sex", "age", "tnm.stage", "tnm.t", "tnm.n", "tnm.m",
      "lymphnodes", "tumor.location", "grade", "chemotherapy.adjuvant",
      "rfs.event", "rfs.delay", "os.event", "os.delay", "dfs.event",
      "dfs.delay", "mmr.status", "cimp.status", "cin.status",
      "tp53.mutation", "kras.mutation", "braf.mutation", "dataset", "cms"
    )
  )
  expect_true(all(
    c("GSE14333", "GSE39582", "TCGA") %in% unique(clin_crc$dataset)
  ))
})
