# load("data-raw/crc_clin_20180713.RData")
# clin_crc <- All.clin
#
# # remove duplicated one
# clin_crc[grep("GSM972295", clin_crc$sample)[2], "sample"] <- "GSM972296"
# colnames(clin_crc) <- tolower(colnames(clin_crc))
#
# # fix cms labels
# clin_crc$cms[clin_crc$cms == "NOLBL"] <- NA
#
# # add rownames for index
# rownames(clin_crc) <- clin_crc$sample
#
# clin_crc[clin_crc$dataset == "GSE14333", c("rfs.event")] <- clin_crc[clin_crc$dataset == "GSE14333", c("dfs.event")]

# load("data/clin_crc.rda")
# write.csv(clin_crc, file = "data-raw/clin_crc.csv")
#
# clin_crc <- read.csv(file = "data-raw/clin_crc.csv", stringsAsFactors = F)[, -1]
# rownames(clin_crc) <- clin_crc$sample
#
# clin_crc$cin.status <- factor(clin_crc$cin.status, levels = c("-", "+"))
#
# devtools::use_data(clin_crc, overwrite = T)


# tcga clinical information
extract_tcga_crc_clinical <- function (clin)
{
  # sample
  sample <- as.character(toupper(clin[, "patient.bcr_patient_barcode"]))
  # sex
  gender <- clin[, "patient.gender"]
  sex <- as.character(factor(gender, levels = c("male", "female"), labels=c("M", "F")))
  # age
  age <- as.numeric(as.character(clin[, "patient.age_at_initial_pathologic_diagnosis"]))
  # tnm.stage
  stage <- clin[, "patient.stage_event.pathologic_stage"]
  stage[grep("iv", stage)] <- 4
  stage[grep("iii", stage)] <- 3
  stage[grep("ii", stage)] <- 2
  stage[grep("i", stage)] <- 1
  tnm.stage <- as.numeric(stage)

  # tnm.t
  Tstage <- clin[, "patient.stage_event.tnm_categories.pathologic_categories.pathologic_t"]
  Tstage[grep("t4", Tstage)] <- 4
  Tstage[grep("t3", Tstage)] <- 3
  Tstage[grep("t2", Tstage)] <- 2
  Tstage[grep("t1", Tstage)] <- 1
  Tstage[grep("tis", Tstage)] <- 1
  tnm.t <- as.numeric(Tstage)

  # tnm.n
  Nstage <- clin[, "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n"]
  Nstage[grep("n2", Nstage)] <- 2
  Nstage[grep("n1", Nstage)] <- 1
  Nstage[grep("n0", Nstage)] <- 0
  tnm.n <- as.numeric(Nstage)

  # tnm.m
  Mstage <- clin[, "patient.stage_event.tnm_categories.pathologic_categories.pathologic_m"]
  Mstage[grep("m1", Mstage)] <- 1
  Mstage[grep("m0", Mstage)] <- 0
  tnm.m <- as.numeric(Mstage)
  # lymphnodes

  # tumor.location
  location <- clin[, "patient.tumor_samples.tumor_sample.tumor_locations.tumor_location.site_of_disease_description"]
  ind <- location %in% c("cecum", "ascending colon", "hepatic flexure", "transverse colon")
  location[ind] <- "R"
  location[!ind] <- "L"
  tumor.location <- location
  # grade

  # chemotherapy.adjuvant
  chemo <- clin[, "patient.drugs.drug.therapy_types.therapy_type"]
  chemotherapy.adjuvant <- rep(0, nrow(clin))
  chemotherapy.adjuvant[!is.na(chemo)] <- 1

  # Survival
  labels <- clin[, "patient.follow_ups.follow_up.new_tumor_events.new_tumor_event_after_initial_treatment"]
  names(labels) <- toupper(clin[, "patient.bcr_patient_barcode"])
  relapse <- (labels == "yes")
  names(relapse) <- names(labels)
  t <- c("patient.follow_ups.follow_up.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment",
         "patient.follow_ups.follow_up.days_to_death", "patient.follow_ups.follow_up.days_to_last_followup")
  tmp <- clin[, t]
  colnames(tmp) <- NULL
  time <- tmp[, 1]
  time[is.na(time)] <- as.vector(t(tmp[is.na(time), 2]))
  time[is.na(time)] <- as.vector(t(tmp[is.na(time), 3]))
  relapse.time <- as.numeric(time)
  names(relapse.time) <- toupper(clin[, "patient.bcr_patient_barcode"])
  labels <- clin[, "patient.vital_status"]
  names(labels) <- toupper(clin[, "patient.bcr_patient_barcode"])
  death <- (labels == "dead")
  names(death) <- names(labels)
  t <- c("patient.follow_ups.follow_up.days_to_death",
         "patient.follow_ups.follow_up.days_to_last_followup")
  tmp <- clin[, t]
  colnames(tmp) <- NULL
  time <- tmp[, 1]
  time[is.na(time)] <- as.vector(t(tmp[is.na(time), 2]))
  death.time <- as.numeric(time)
  names(death.time) <- toupper(clin[, "patient.bcr_patient_barcode"])
  DFS <- relapse | death
  DFS[is.na(DFS) & death] <- T

  # rfs.event
  rfs.event <- as.integer(relapse)
  # rfs.delay
  rfs.delay <- relapse.time/30
  # os.event
  os.event <- as.integer(death)
  # os.delay
  os.delay <- death.time/30
  # dfs.event
  dfs.event <- as.integer(DFS)
  # dfs.delay
  dfs.delay <- rfs.delay
  # mmr.status
  mmr.raw <- clin[, "patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status"]
  mmr <- rep(NA, nrow(clin))
  mmr[mmr.raw == "mss"] <- "mss"
  mmr[grep("msi", mmr.raw)] <- "msi"
  mmr.status <- mmr

  # cimp.status
  # cin.status
  # tp53.mutation
  # kras.mutation
  kras.mutation <- as.integer(clin[, "patient.kras_mutation_found"] == "yes")
  # braf.mutation
  braf.mutation <- as.integer(clin[, "patient.braf_gene_analysis_result"] == "abnormal")
  # dataset
  # cms

  data.frame(sample, sex, age, tnm.stage, tnm.t, tnm.n, tnm.m, lymphnodes=NA, tumor.location, grade=NA, chemotherapy.adjuvant,
             rfs.event, rfs.delay, os.event, os.delay, dfs.event, dfs.delay,
             mmr.status, cimp.status=NA, cin.status=NA, tp53.mutation=NA, kras.mutation, braf.mutation, dataset="TCGA")

}

# library(RTCGA.clinical.20160128)
# tcga.clin.raw <- COADREAD.clinical.20160128
# clin.tcga <- extract_tcga_crc_clinical(tcga.clin.raw)
#
# # relapce CMS etc.
# load("data-raw/clin_crc_gf.RData")
# clin.tcga["TCGA-D5-5537", "tnm.stage"] <- 3
# clin.tcga$cms <- clin_crc_gf[rownames(clin.tcga), "cms"]
#
# clin_crc_gf <- rbind(clin.tcga, clin_crc_gf[clin_crc_gf$dataset != "TCGA", ])
#
# devtools::use_data(clin_crc_gf, overwrite = T)
#
# clin_crc_gf[clin_crc_gf$dataset == "GSE14333", "rfs.event"] <- 1-clin_crc_gf[clin_crc_gf$dataset == "GSE14333", "rfs.event"]
# clin_crc_gf[clin_crc_gf$dataset == "GSE14333", "dfs.event"] <- 1-clin_crc_gf[clin_crc_gf$dataset == "GSE14333", "dfs.event"]
#
# devtools::use_data(clin_crc_gf, overwrite = T)

##-----------------------------------------
## 2019.6.10


# load("data/clin_crc_gf.rda")
# write.csv(clin_crc_gf, file = "data-raw/clin_crc_gf.csv")

clin_crc_gf <- read.csv(file = "data-raw/clin_crc_gf.csv", stringsAsFactors = F)[, -1]
rownames(clin_crc_gf) <- clin_crc_gf$sample

tcga_cdr_crc <- read.csv(file = "data-raw/tcga_cdr_crc.csv", stringsAsFactors = F)
rownames(tcga_cdr_crc) <- tcga_cdr_crc[, 1]

pt_ids <- intersect(clin_crc_gf$sample, rownames(tcga_cdr_crc))
clin_crc_gf[pt_ids, "dfs.delay"] <- as.numeric(tcga_cdr_crc[pt_ids, "PFS.time"])/30
clin_crc_gf[pt_ids, "dfs.event"] <- as.numeric(tcga_cdr_crc[pt_ids, "PFS"])

clin_crc_gf[pt_ids, "rfs.delay"] <- as.numeric(tcga_cdr_crc[pt_ids, "PFI.time.1"])/30
clin_crc_gf[pt_ids, "rfs.event"] <- as.numeric(tcga_cdr_crc[pt_ids, "PFI.1"])

clin_crc_gf[pt_ids, "os.delay"] <- as.numeric(tcga_cdr_crc[pt_ids, "DSS.time.cr"])/30
clin_crc_gf[pt_ids, "os.event"] <- as.numeric(as.numeric(tcga_cdr_crc[pt_ids, "DSS_cr"])>=1)

clin_crc_cell <- clin_crc_gf

devtools::use_data(clin_crc_cell, overwrite = T)
