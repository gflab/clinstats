toy_tcga_table <- function() {
  data.frame(
    patient.bcr_patient_barcode = c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"),
    patient.age_at_initial_pathologic_diagnosis = c(60, 70, 55),
    patient.gender = c("male", "female", "female"),
    patient.vital_status = c("dead", "alive", "alive"),
    patient.stage_event.pathologic_stage = c("stage iiib", "stage i", "stage iv"),
    patient.stage_event.tnm_categories.pathologic_categories.pathologic_t = c("t3", "t1", "t4"),
    patient.stage_event.tnm_categories.pathologic_categories.pathologic_n = c("n1", "n0", "nx"),
    patient.follow_ups.follow_up.primary_therapy_outcome_success = c(
      "partial remission/response", NA, "stable disease"
    ),
    patient.drugs.drug.therapy_types.therapy_type = c("chemotherapy", NA, "chemotherapy"),
    patient.follow_ups.follow_up.new_tumor_events.new_tumor_event_after_initial_treatment = c(
      "yes", NA, NA
    ),
    patient.follow_ups.follow_up.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment = c(
      200, NA, NA
    ),
    patient.follow_ups.follow_up.days_to_death = c(500, NA, NA),
    patient.follow_ups.follow_up.days_to_last_followup = c(NA, 300, 800),
    check.names = FALSE
  )
}

test_that("tcga_clinical maps the default layout", {
  result <- tcga_clinical(toy_tcga_table())

  expect_s3_class(result, "tbl_df")
  expect_equal(result$sample, c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"))
  expect_equal(result$stage, c(3, 1, 4))
  expect_equal(result$t_stage, c(3, 1, 4))
  expect_equal(result$node_positive, c(TRUE, FALSE, NA))
  expect_equal(result$relapse, c(TRUE, NA, NA))
  expect_equal(result$relapse_time, c(200, 300, 800))
  expect_equal(result$death, c(TRUE, FALSE, FALSE))
  expect_equal(result$death_time, c(500, 300, 800))
  expect_equal(result$dfs, c(TRUE, NA, NA))
  expect_equal(result$chemotherapy, c(1L, 0L, 1L))
  expect_equal(result$chemo_response, c("p_sen", NA, "resist"))
})

test_that("stage mapping distinguishes roman numerals", {
  values <- c("stage i", "stage ia", "stage ii", "stage iia", "stage iiib", "stage iv", "III")
  expect_equal(map_stage(values), c(1, 1, 2, 2, 3, 4, 3))
  expect_equal(map_stage(c(NA, "")), c(NA_real_, NA_real_))
})

test_that("tcga_clinical reports missing columns", {
  table <- toy_tcga_table()
  table[["patient.vital_status"]] <- NULL
  expect_error(tcga_clinical(table), "vital_status")
})
