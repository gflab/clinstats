#' Extract clinical variables from a legacy TCGA clinical table
#'
#' Converts the wide clinical export distributed with legacy TCGA data
#' releases (column names of the form
#' `patient.age_at_initial_pathologic_diagnosis`) into an analysis-ready
#' tibble with follow-up times in days. Three column layouts are supported:
#' the generic layout, the oesophageal carcinoma layout (`"esca"`), and the
#' stomach adenocarcinoma layout (`"stad"`).
#'
#' @param clin Data frame with one row per patient in the legacy TCGA layout.
#' @param cancer Column layout: `"default"`, `"esca"` (the historical
#'   `ESAC = TRUE` variant), or `"stad"`.
#'
#' @return A tibble with the patient barcode, age, sex, stage, node status,
#'   T stage, disease-free survival, relapse and death status with their
#'   times, and, where available, chemotherapy and response. The `"esca"`
#'   layout additionally returns the histological subtype and tumour
#'   location.
#'
#' @export
tcga_clinical <- function(clin, cancer = c("default", "esca", "stad")) {
  cancer <- match.arg(cancer)
  if (!is.data.frame(clin)) {
    cli::cli_abort("{.arg clin} must be a data frame.")
  }

  common <- c(
    "patient.bcr_patient_barcode",
    "patient.gender",
    "patient.vital_status",
    "patient.stage_event.pathologic_stage",
    "patient.stage_event.tnm_categories.pathologic_categories.pathologic_t",
    "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n"
  )
  follow_up <- "patient.follow_ups.follow_up"
  required <- switch(cancer,
    default = c(
      common,
      "patient.age_at_initial_pathologic_diagnosis",
      paste0(follow_up, ".primary_therapy_outcome_success"),
      "patient.drugs.drug.therapy_types.therapy_type",
      paste0(follow_up, ".new_tumor_events.new_tumor_event_after_initial_treatment"),
      paste0(follow_up, ".new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment"),
      paste0(follow_up, ".days_to_death"),
      paste0(follow_up, ".days_to_last_followup")
    ),
    esca = c(
      common,
      "patient.primary_pathology.age_at_initial_pathologic_diagnosis",
      "patient.new_tumor_events.new_tumor_event_after_initial_treatment",
      "patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment",
      "patient.days_to_death",
      paste0(follow_up, ".days_to_last_followup"),
      "patient.biospecimen_cqcf.histological_type",
      "patient.primary_pathology.esophageal_tumor_cental_location"
    ),
    stad = c(
      common,
      "patient.age_at_initial_pathologic_diagnosis",
      paste0(follow_up, ".new_tumor_event_after_initial_treatment"),
      paste0(follow_up, ".days_to_new_tumor_event_after_initial_treatment"),
      paste0(follow_up, ".days_to_death"),
      paste0(follow_up, ".days_to_last_followup")
    )
  )
  check_columns(clin, required)

  col <- function(name) clin[[name]]
  sample <- toupper(as.character(col("patient.bcr_patient_barcode")))
  age <- as.numeric(if (cancer == "esca") {
    col("patient.primary_pathology.age_at_initial_pathologic_diagnosis")
  } else {
    col("patient.age_at_initial_pathologic_diagnosis")
  })
  sex <- as.character(col("patient.gender"))
  stage <- map_stage(col("patient.stage_event.pathologic_stage"))
  t_stage <- map_t_stage(
    col("patient.stage_event.tnm_categories.pathologic_categories.pathologic_t")
  )
  n_stage <- tolower(as.character(
    col("patient.stage_event.tnm_categories.pathologic_categories.pathologic_n")
  ))
  node_positive <- !grepl("^n0", n_stage)
  node_positive[grepl("^nx", n_stage) | is.na(n_stage)] <- NA

  relapse_label <- switch(cancer,
    default = col(paste0(
      follow_up,
      ".new_tumor_events.new_tumor_event_after_initial_treatment"
    )),
    esca = col("patient.new_tumor_events.new_tumor_event_after_initial_treatment"),
    stad = col(paste0(follow_up, ".new_tumor_event_after_initial_treatment"))
  )
  relapse <- as.character(relapse_label) == "yes"
  relapse[is.na(relapse_label)] <- NA

  relapse_time <- switch(cancer,
    default = coalesce_numeric(
      col(paste0(
        follow_up,
        ".new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment"
      )),
      col(paste0(follow_up, ".days_to_death")),
      col(paste0(follow_up, ".days_to_last_followup"))
    ),
    esca = coalesce_numeric(
      col("patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment"),
      col("patient.days_to_death"),
      col(paste0(follow_up, ".days_to_last_followup"))
    ),
    stad = coalesce_numeric(
      col(paste0(follow_up, ".days_to_new_tumor_event_after_initial_treatment")),
      col(paste0(follow_up, ".days_to_death")),
      col(paste0(follow_up, ".days_to_last_followup"))
    )
  )

  vital <- tolower(as.character(col("patient.vital_status")))
  death <- vital == "dead"
  death[is.na(vital)] <- NA
  death_time <- coalesce_numeric(
    if (cancer == "esca") {
      col("patient.days_to_death")
    } else {
      col(paste0(follow_up, ".days_to_death"))
    },
    col(paste0(follow_up, ".days_to_last_followup"))
  )
  dfs <- relapse | death
  dfs[is.na(dfs) & !is.na(death) & death] <- TRUE

  out <- tibble::tibble(
    sample = sample,
    age = age,
    sex = sex,
    stage = stage,
    node_positive = node_positive,
    t_stage = t_stage,
    dfs = dfs,
    relapse = relapse,
    relapse_time = relapse_time,
    death = death,
    death_time = death_time
  )

  if (cancer == "default") {
    chemo <- as.character(col("patient.drugs.drug.therapy_types.therapy_type"))
    response <- as.character(col(
      paste0(follow_up, ".primary_therapy_outcome_success")
    ))
    out$chemotherapy <- ifelse(is.na(chemo), 0, as.integer(chemo == "chemotherapy"))
    out$chemo_response <- map_chemo_response(response)
  }
  if (cancer == "esca") {
    histology <- tolower(as.character(col("patient.biospecimen_cqcf.histological_type")))
    out$subtype <- ifelse(
      histology == "esophagus adenocarcinoma  nos",
      "EAC",
      "ESCC"
    )
    out$subtype[is.na(histology)] <- NA
    out$location <- as.character(col(
      "patient.primary_pathology.esophageal_tumor_cental_location"
    ))
  }
  out
}

map_stage <- function(value) {
  text <- tolower(trimws(as.character(value)))
  roman <- text
  has_stage <- grepl("stage", text)
  roman[has_stage] <- sub("^.*stage\\s*([ivx]+).*$", "\\1", text[has_stage])
  roman[!has_stage] <- gsub("[^ivx]", "", roman[!has_stage])
  out <- rep(NA_real_, length(text))
  out[roman == "i"] <- 1
  out[roman == "ii"] <- 2
  out[roman == "iii"] <- 3
  out[roman == "iv"] <- 4
  out
}

map_t_stage <- function(value) {
  text <- tolower(as.character(value))
  out <- rep(NA_real_, length(text))
  out[grepl("^tis", text)] <- 1
  out[grepl("^t1", text)] <- 1
  out[grepl("^t2", text)] <- 2
  out[grepl("^t3", text)] <- 3
  out[grepl("^t4", text)] <- 4
  out
}

map_chemo_response <- function(value) {
  text <- tolower(trimws(value))
  out <- rep(NA_character_, length(text))
  out[text == "complete remission/response"] <- "sen"
  out[text == "partial remission/response"] <- "p_sen"
  out[text %in% c("stable disease", "progressive disease")] <- "resist"
  out
}

coalesce_numeric <- function(...) {
  columns <- list(...)
  out <- as.numeric(columns[[1]])
  if (length(columns) > 1) {
    for (x in columns[-1]) {
      x <- as.numeric(x)
      out[is.na(out)] <- x[is.na(out)]
    }
  }
  out
}
