# ADAPTED FILE 01: advanced CKD / KRT versus eGFR >60 at the case index date.
# This adapted script does not overwrite the old all-stage CKD database caches.
# eGFR source: all_patid / clean_egfr_medcodes (patid, date, testvalue),
# identified by supplied file 02. The clean cache must already exist.
# ALL numeric results are retained until selecting the latest pre-index result.
# The six-monthly baseline summaries in 02 are not used for case-specific dates.
#
# Expensive query outputs are cached BEFORE collect(). Local matching tables
# remain available; final cohorts are also cached and exposed as *_db objects.
# Local codelist joins use temporary database copies as in file 04.
# Completed persistent caches are reusable after reconnecting. This does not
# keep an R process alive or guarantee completion of an interrupted SQL write.
#
# Design choices retained / made explicit:
# - Cases: earliest recorded stage_4, stage_5, or qualifying local code date.
# - Stage dates are trusted as supplied: inspect the upstream algorithm to
#   confirm chronicity and whether dates refer to first evidence or confirmation.
# - The user confirms that the dialysis codelists identify established KRT,
#   not acute treatment or preparation for KRT.
# - Controls: never recorded as ADVANCED CKD / KRT in the available data;
#   this retains your earlier decision to exclude future cases, but permits
#   mild CKD. It uses future exposure information and is a study-design choice.
# - Latest pre-index eGFR >60 within 365 days; no qualifying result = ineligible.
#   This is a recent preserved-eGFR comparator, not proof of no CKD or of
#   persistently preserved eGFR. It can include recovered lower eGFR.
# - Age >18; exact practice; nearest DOB, gender, registration start/end;
#   aim for four controls, retain one to four, without replacement.
# - Six months BEFORE and AFTER index, as actually implemented in uploaded 01.
# - Membership in diabetes_cohort does not establish diabetes before index.
# - Exclude cases and candidate controls with any valid alldementia diagnosis
#   on or before the CASE index date + six calendar months (inclusive).
#   Later dementia is retained as a potential outcome, not an exclusion.
# - This defines a cohort dementia-free through six months. Outcome follow-up
#   and person-time should be defined consistently with that eligibility period;
#   this script retains the original kidney index_date and does not model outcomes.
# - A post-index observation requirement conditions on future follow-up.

############################################################################################

#Setup
library(dplyr)
library(tidyr)
library(lubridate)
library(MatchIt)
library(tidyverse)
library(aurum)
library(EHRBiomarkr)


cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes_2024 = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix = "ckd"

# Load Robert's additional codelists locally

codelist_root <- paste0(
  "C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/",
  "Github clone/CPRDKidneyFailureCognitive/CPRD-Codelists"
)

custom_codelist_directories <- c(
  file.path(codelist_root, "Medcodes"),
  file.path(codelist_root, "ICD10"),
  file.path(codelist_root, "OPCS4")
)

missing_directories <- custom_codelist_directories[
  !dir.exists(custom_codelist_directories)
]

if (length(missing_directories) > 0) {
  stop(
    "The following codelist directories were not found: ",
    paste(missing_directories, collapse = ", ")
  )
}

read_local_codelists <- function(directories) {

  files <- unlist(
    lapply(
      directories,
      list.files,
      pattern = "\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
  )

  output <- list()

  for (file in files) {

    codelist <- readr::read_tsv(
      file,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
      rename_with(stringr::str_to_lower)

    code_columns <- intersect(
      names(codelist),
      c("medcodeid", "icd10", "opcs4")
    )

    if (length(code_columns) != 1) {
      warning(
        "Skipping ", file,
        ": expected exactly one of medcodeid, icd10 or opcs4"
      )
      next
    }

    code_column <- code_columns[[1]]

    codelist <- codelist %>%
      filter(
        !is.na(.data[[code_column]]),
        .data[[code_column]] != ""
      )

    code_name <- file %>%
      basename() %>%
      tools::file_path_sans_ext() %>%
      stringr::str_to_lower() %>%
      stringr::str_remove("^exeter_medcodelist_") %>%
      stringr::str_remove("^exeter_")

    if (code_name %in% names(output)) {
      stop(
        "More than one local codelist generated the name: ",
        code_name
      )
    }

    output[[code_name]] <- codelist
  }

  output
}

custom_codes <- read_local_codelists(
  custom_codelist_directories
)

sort(names(custom_codes))


# Settings: review these before running.
pre_index_months <- 6L
post_index_months <- 6L
dementia_lag_months <- 6L
egfr_lookback_days <- 365L       # Proposed window; change to your chosen window.
hes_end_date <- as.Date("2023-03-31") # Same linked-data cutoff as supplied 04.

# Cache naming: keep this version unchanged to REUSE completed tables.
# Change it when changing source data, codelists, validity rules or cohort /
# matching settings. In particular, final cached outputs are NOT overwritten
# merely because you rerun matching with different settings or another seed.
cache_version <- "v1"
cache_name <- function(stem) paste0("rk_adv_", cache_version, "_", stem)
input_analysis <- cprd$analysis("all_patid")
output_analysis <- cprd$analysis("ckd")
database_source <- dbplyr::remote_src(cprd$tables$patient)

# For already-local results only. Expensive SQL inputs are cached directly below.
cache_local_output <- function(data, stem, indexes) {
  id_columns <- intersect(c("patid", "ckd_patid", "pracid"), names(data))
  id_types <- setNames(rep("VARCHAR(64)", length(id_columns)), id_columns)
  uploaded <- dplyr::copy_to(
    database_source, data,
    name = paste0("tmp_", cache_name(stem)),
    types = id_types, temporary = TRUE, overwrite = TRUE,
    analyze = FALSE, in_transaction = FALSE
  )
  uploaded %>% output_analysis$cached(
    cache_name(stem), unique_indexes = "patid", indexes = indexes
  )
}

# 1. Earliest algorithm stage 4 and stage 5 dates.
analysis <- cprd$analysis("all_patid")
ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  analysis$cached("ckd_stages_from_algorithm")

algorithm_dates <- ckd_stages_from_algorithm %>%
  select(patid, stage_4, stage_5) %>%
  filter(!is.na(stage_4) | !is.na(stage_5)) %>%
  input_analysis$cached(cache_name("algorithm_stage45"), indexes = "patid") %>%
  collect() %>%
  mutate(patid = as.character(patid)) %>%
  pivot_longer(c(stage_4, stage_5), names_to = "criterion", values_to = "date") %>%
  filter(!is.na(date)) %>%
  mutate(date = as.Date(date), source = "egfr_algorithm") %>%
  group_by(patid, criterion, source) %>%
  summarise(date = min(date), .groups = "drop")

# 2. Extract five local advanced CKD / KRT conditions plus alldementia, following 04.
# The loader lowercases filenames, so ckd5_noKRT becomes ckd5_nokrt.
advanced_conditions <- c("ckd5_nokrt", "ckd5", "haemodialysis",
                         "peritoneal_dialysis", "transplant")

first_valid_date <- function(events, source_name, stem) {
  events %>%
    inner_join(cprd$tables$validDateLookup, by = "patid") %>%
    filter(!is.na(date), date >= min_dob) %>%
    { if (source_name == "gp") filter(., date <= gp_end_date)
      else filter(., date <= !!hes_end_date) } %>%
    group_by(patid) %>%
    summarise(date = min(date, na.rm = TRUE), .groups = "drop") %>%
    input_analysis$cached(
      cache_name(paste0("first_", stem)),
      unique_indexes = "patid", indexes = "date"
    ) %>%
    collect() %>%
    mutate(patid = as.character(patid), date = as.Date(date), source = source_name)
}

code_dates_list <- list()
for (condition in c(advanced_conditions, "alldementia")) {
  found <- FALSE
  med <- custom_codes[[condition]]
  icd <- custom_codes[[paste0("icd10_", condition)]]
  opc <- custom_codes[[paste0("opcs4_", condition)]]

  if (!is.null(med) && nrow(med) > 0L) {
    found <- TRUE
    events <- cprd$tables$observation %>%
      inner_join(distinct(med, medcodeid), by = "medcodeid", copy = TRUE) %>%
      select(patid, date = obsdate) %>%
      input_analysis$cached(
        cache_name(paste0("raw_", condition, "_gp")),
        indexes = c("patid", "date")
      )
    code_dates_list[[paste0(condition, "_gp")]] <-
      first_valid_date(events, "gp", paste0(condition, "_gp")) %>%
      mutate(criterion = condition)
  }
  if (!is.null(icd) && nrow(icd) > 0L) {
    found <- TRUE
    events <- cprd$tables$hesDiagnosisEpi %>%
      inner_join(distinct(icd, icd10),
                 sql_on = "LHS.ICD LIKE CONCAT(icd10,'%')", copy = TRUE) %>%
      select(patid, date = epistart) %>%
      input_analysis$cached(
        cache_name(paste0("raw_", condition, "_icd10")),
        indexes = c("patid", "date")
      )
    code_dates_list[[paste0(condition, "_icd10")]] <-
      first_valid_date(events, "hes", paste0(condition, "_icd10")) %>%
      mutate(criterion = condition)
  }
  if (!is.null(opc) && nrow(opc) > 0L) {
    found <- TRUE
    events <- cprd$tables$hesProceduresEpi %>%
      inner_join(distinct(opc, opcs4), by = c("OPCS" = "opcs4"), copy = TRUE) %>%
      select(patid, date = evdate) %>%
      input_analysis$cached(
        cache_name(paste0("raw_", condition, "_opcs4")),
        indexes = c("patid", "date")
      )
    code_dates_list[[paste0(condition, "_opcs4")]] <-
      first_valid_date(events, "hes", paste0(condition, "_opcs4")) %>%
      mutate(criterion = condition)
  }
  if (!found) stop("No nonempty local codelist found for: ", condition)
}

coded_evidence <- bind_rows(code_dates_list)

# Keep dementia evidence separate: it must NEVER contribute to a kidney index.
dementia_dates <- coded_evidence %>%
  filter(criterion == "alldementia") %>%
  group_by(patid) %>%
  summarise(first_dementia_date = min(date), .groups = "drop")

advanced_evidence <- bind_rows(
  algorithm_dates,
  coded_evidence %>% filter(criterion %in% advanced_conditions)
) %>%
  distinct(patid, date, criterion, source)

# Audit overlapping sources; these counts are NOT mutually exclusive.
advanced_evidence %>% count(criterion, source) %>% print(n = Inf, width = Inf)

advanced_ckd_ids <- advanced_evidence %>%
  group_by(patid) %>%
  summarise(first_advanced_ckd_date = min(date), .groups = "drop")

# Keep evidence supporting the entry date. Multiple criteria can tie.
advanced_index_evidence <- advanced_evidence %>%
  inner_join(advanced_ckd_ids, by = "patid") %>%
  filter(date == first_advanced_ckd_date)

# Save derived ID/date tables before the long matching stage.
advanced_ckd_ids_db <- cache_local_output(
  advanced_ckd_ids, "advanced_ckd_ids", indexes = "first_advanced_ckd_date"
)
dementia_dates_db <- cache_local_output(
  dementia_dates, "dementia_dates", indexes = "first_dementia_date"
)

# 3. Apply the existing practice/gender exclusions to all diabetes patients.
analysis <- cprd$analysis("diabetes_cohort")
practice_exclusion_ids <- practice_exclusion_ids %>%
  analysis$cached("practice_exclusion_ids")
gender_exclusion_ids <- gender_exclusion_ids %>%
  analysis$cached("gender_exclusion_ids")
analysis <- cprd$analysis("all")
diabetes_cohort <- diabetes_cohort %>% analysis$cached("diabetes_cohort")

matching_columns <- c("patid", "pracid", "dob", "gender",
                      "regstartdate", "gp_end_date")
diabetes_matching <- diabetes_cohort %>%
  anti_join(practice_exclusion_ids, by = "patid") %>%
  anti_join(gender_exclusion_ids, by = "patid") %>%
  select(all_of(matching_columns)) %>%
  input_analysis$cached(cache_name("diabetes_matching_base"),
                        unique_indexes = "patid", indexes = "pracid") %>%
  collect() %>%
  mutate(patid = as.character(patid)) %>%
  left_join(dementia_dates, by = "patid")

rk_diabetes_advanced_ckd_cohort <- diabetes_matching %>%
  inner_join(advanced_ckd_ids, by = "patid")

rk_diabetes_advanced_ckd_cohort_db <- cache_local_output(
  rk_diabetes_advanced_ckd_cohort, "diabetes_advanced_ckd_cohort",
  indexes = c("pracid", "first_advanced_ckd_date")
)

cases_raw <- rk_diabetes_advanced_ckd_cohort %>%
  rename(index_date = first_advanced_ckd_date)

# Deliberately exclude advanced_ckd_ids, NOT the old all-stage ckd_ids.
controls_raw <- diabetes_matching %>%
  anti_join(advanced_ckd_ids, by = "patid")

print(tibble(advanced_ckd_cases = nrow(cases_raw),
             possible_comparators_before_egfr_check = nrow(controls_raw)),
      width = Inf)

# 4. Load the clean dated eGFR source identified in file 02.
# Follow the existing aurum cached-table loading convention.
analysis <- cprd$analysis("all_patid")
clean_egfr_medcodes <- clean_egfr_medcodes %>%
  analysis$cached("clean_egfr_medcodes")

egfr_history <- clean_egfr_medcodes %>%
  select(patid, egfr_date = date, egfr = testvalue) %>%
  filter(!is.na(egfr_date), !is.na(egfr), egfr > 0) %>%
  group_by(patid, egfr_date) %>%
  summarise(egfr = min(egfr, na.rm = TRUE), .groups = "drop") %>%
  input_analysis$cached(cache_name("egfr_daily"),
                        indexes = c("patid", "egfr_date")) %>%
  collect()
stopifnot(is.numeric(egfr_history$egfr))

# Use latest pre-index measurement within 365 days, preserving the rule in
# the previous matching draft. This differs intentionally from file 02's
# nearest -730/+7-day baseline summaries. We also retain the conservative
# same-day minimum rule instead of file 02's maximum for eGFR.
egfr_data <- egfr_history %>%
  mutate(patid = as.character(patid), egfr_date = as.Date(egfr_date)) %>%
  semi_join(controls_raw %>% select(patid), by = "patid") %>%
  filter(!is.na(egfr_date), is.finite(egfr), egfr > 0) %>%
  # Same-day minimum has already been calculated and cached on the server.
  arrange(patid, egfr_date)

# Reduce the candidate pool, but DO NOT filter the history to only high values.
controls_raw <- controls_raw %>%
  semi_join(egfr_data %>% filter(egfr > 60) %>% distinct(patid), by = "patid")

egfr_rows <- split(seq_len(nrow(egfr_data)), egfr_data$patid)

latest_egfr <- function(patient_id, index_dates) {
  rows <- egfr_rows[[as.character(patient_id)]]
  values <- rep(NA_real_, length(index_dates))
  dates <- as.Date(rep(NA_character_, length(index_dates)))
  if (length(rows)) {
    h <- egfr_data[rows, , drop = FALSE]
    pos <- findInterval(as.numeric(index_dates), as.numeric(h$egfr_date))
    found <- pos > 0L
    values[found] <- h$egfr[pos[found]]
    dates[found] <- h$egfr_date[pos[found]]
  }
  tibble(egfr_at_index = values, egfr_date = dates)
}

# Do not silently discard repeated patients or conflicting demographics.
stopifnot(!anyDuplicated(cases_raw$patid),
          !anyDuplicated(controls_raw$patid))

prepare_patients <- function(x) {
  x %>%
    mutate(patid = as.character(patid), pracid = as.character(pracid),
           gender = as.character(gender),
           across(any_of(c("dob", "regstartdate", "gp_end_date", "index_date",
                             "first_dementia_date")),
                  as.Date)) %>%
    drop_na(all_of(c("patid", "pracid", "dob", "gender",
                     "regstartdate", "gp_end_date"))) %>%
    filter(regstartdate <= gp_end_date) %>%
    arrange(patid)
}

cases <- prepare_patients(cases_raw) %>%
  drop_na(index_date) %>%
  mutate(
    window_start = index_date %m-% months(pre_index_months),
    window_end = index_date %m+% months(post_index_months),
    dementia_exclusion_end = index_date %m+% months(dementia_lag_months)
  ) %>%
  filter(
    is.na(first_dementia_date) | first_dementia_date > dementia_exclusion_end,
    index_date > dob %m+% years(18),
    regstartdate <= window_start,
    gp_end_date >= window_end
  )

controls <- prepare_patients(controls_raw)
stopifnot(!any(cases$patid %in% controls$patid))

# Empty result with stable column types, including when no matches are found.
empty_pairs <- tibble(ckd_patid = character(), patid = character(),
                      index_date = as.Date(character()), pracid = character())

# Optional maximum age difference:
# Inf = no limit; change to 5 for a maximum gap of approximately five years.
max_age_gap_years <- Inf

match_one_practice <- function(ca, co) {

  # Only one available control is now needed.
  if (nrow(ca) == 0L || nrow(co) == 0L) return(empty_pairs)

  # Each control must span the CASE'S registration window
  # and be older than 18 at that case's index date.
  eligible <- outer(
    as.numeric(ca$window_start),
    as.numeric(co$regstartdate),
    `>=`
  ) &
    outer(
      as.numeric(ca$window_end),
      as.numeric(co$gp_end_date),
      `<=`
    ) &
    outer(
      as.numeric(ca$index_date),
      as.numeric(co$dob %m+% years(18)),
      `>`
    )

  # Dementia eligibility also depends on the CASE'S index date.
  # No recorded diagnosis gets Inf so it remains eligible. Any first diagnosis
  # on/before index + six months disallows that pair; later diagnoses are allowed.
  control_dementia_day <- as.numeric(co$first_dementia_date)
  control_dementia_day[is.na(control_dementia_day)] <- Inf
  eligible <- eligible & outer(
    as.numeric(ca$dementia_exclusion_end), control_dementia_day, `<`
  )

  # Each control needs a RECENT latest eGFR >60 at THIS case's index.
  # A historical high result is insufficient if a newer pre-index result is low.
  for (j in seq_len(nrow(co))) {
    lab <- latest_egfr(co$patid[j], ca$index_date)
    recent <- !is.na(lab$egfr_at_index) & !is.na(lab$egfr_date) &
      lab$egfr_at_index > 60 &
      lab$egfr_date >= ca$index_date - egfr_lookback_days
    eligible[, j] <- eligible[, j] & recent
  }

  # Optional hard limit on case-control age differences.
  if (is.finite(max_age_gap_years)) {
    age_gap <- abs(
      outer(as.numeric(ca$dob), as.numeric(co$dob), `-`)
    ) / 365.25

    eligible <- eligible & age_gap <= max_age_gap_years
  }

  # Retain cases with at least ONE potentially eligible control.
  keep_cases <- rowSums(eligible) >= 1L
  ca <- ca[keep_cases, , drop = FALSE]
  eligible <- eligible[keep_cases, , drop = FALSE]

  if (nrow(ca) == 0L) return(empty_pairs)

  keep_controls <- colSums(eligible) > 0L
  co <- co[keep_controls, , drop = FALSE]
  eligible <- eligible[, keep_controls, drop = FALSE]

  dat <- bind_rows(
    mutate(ca, CKD = 1L),
    mutate(co, CKD = 0L)
  ) %>%
    mutate(
      gender = factor(gender),
      dob_number = as.numeric(dob),
      regstart_number = as.numeric(regstartdate),
      gpend_number = as.numeric(gp_end_date)
    ) %>%
    as.data.frame()

  rownames(dat) <- dat$patid

  vars <- c(
    "dob_number", "gender",
    "regstart_number", "gpend_number"
  )

  vars <- vars[
    vapply(dat[vars], function(x) n_distinct(x) > 1L, logical(1))
  ]

  form <- if (length(vars)) {
    reformulate(vars, response = "CKD")
  } else {
    CKD ~ 1
  }

  # With just one case and one control there is no ranking choice;
  # avoid estimating a pooled covariance from two single-person groups.
  if (nrow(ca) == 1L && nrow(co) == 1L) {
    distances <- matrix(0, 1L, 1L)
  } else if (length(vars)) {
    distances <- MatchIt::mahalanobis_dist(form, data = dat)
  } else {
    distances <- matrix(0, nrow(ca), nrow(co))
  }

  dimnames(distances) <- list(ca$patid, co$patid)
  stopifnot(all(is.finite(distances)))

  # Forbid pairs that fail registration, adulthood or optional age-gap rules.
  distances[!eligible] <- Inf

  fit <- matchit(
    form,
    data = dat,
    method = "nearest",
    distance = distances,
    exact = ~ pracid,
    ratio = min(4L, nrow(co)),
    replace = FALSE,
    m.order = "random"
  )

  # Retain every case with at least one actual match.
  mm <- fit$match.matrix
  mm <- mm[rowSums(!is.na(mm)) >= 1L, , drop = FALSE]

  if (nrow(mm) == 0L) return(empty_pairs)

  # Extract the available matches and remove unfilled slots.
  tibble(
    ckd_patid = rep(rownames(mm), each = ncol(mm)),
    patid = as.vector(t(mm))
  ) %>%
    filter(!is.na(patid)) %>%
    left_join(
      ca %>% select(ckd_patid = patid, index_date, pracid),
      by = "ckd_patid"
    )
}

# Work within practices to avoid an enormous whole-cohort distance matrix.
# Large individual practices can still need substantial memory.
case_groups <- split(seq_len(nrow(cases)), cases$pracid)
control_groups <- split(seq_len(nrow(controls)), controls$pracid)
pair_list <- vector("list", length(case_groups))
set.seed(123)

for (i in seq_along(case_groups)) {
  practice <- names(case_groups)[i]
  ca <- cases[case_groups[[i]], , drop = FALSE]
  co_rows <- control_groups[[practice]]
  co <- controls[if (is.null(co_rows)) integer() else co_rows, , drop = FALSE]
  pair_list[[i]] <- match_one_practice(ca, co)
  if (i %% 100L == 0L) message("Processed ", i, " practices")
}

matched_pairs <- bind_rows(empty_pairs, bind_rows(pair_list)) %>%
  add_count(ckd_patid, name = "n_controls")

matched_pairs %>%
  distinct(ckd_patid, n_controls) %>%
  count(n_controls, name = "number_of_ckd_cases")

rk_diabetes_egfr60_matched <- matched_pairs %>%
  left_join(controls %>% select(-pracid), by = "patid")

# Attach the actual result supporting each selected control's eligibility.
# Each control is used only once, so there is one lookup per output row.
control_labs <- lapply(seq_len(nrow(rk_diabetes_egfr60_matched)), function(i) {
  latest_egfr(rk_diabetes_egfr60_matched$patid[i],
              rk_diabetes_egfr60_matched$index_date[i])
})
control_labs <- bind_rows(
  tibble(egfr_at_index = double(), egfr_date = as.Date(character())),
  bind_rows(control_labs)
)
rk_diabetes_egfr60_matched <- bind_cols(rk_diabetes_egfr60_matched, control_labs)
stopifnot(all(rk_diabetes_egfr60_matched$egfr_at_index > 60),
          all(rk_diabetes_egfr60_matched$egfr_date <=
                rk_diabetes_egfr60_matched$index_date),
          all(rk_diabetes_egfr60_matched$egfr_date >=
                rk_diabetes_egfr60_matched$index_date - egfr_lookback_days))

rk_diabetes_advanced_ckd_matched <- cases %>%
  semi_join(matched_pairs, by = c("patid" = "ckd_patid"))

matching_summary <- tibble(
  original_ckd_patients = nrow(cases_raw),
  ckd_patients_excluded_by_dementia_rule = sum(
    !is.na(cases_raw$first_dementia_date) &
      cases_raw$first_dementia_date <=
        as.Date(cases_raw$index_date) %m+% months(dementia_lag_months),
    na.rm = TRUE
  ),
  eligible_ckd_patients_after_all_filters = nrow(cases),
  matched_ckd_patients = nrow(rk_diabetes_advanced_ckd_matched),
  matched_controls = nrow(rk_diabetes_egfr60_matched)
)
print(matching_summary, width = Inf)

# Structural checks: four unique controls, exact practice, shared case date,
# and control registration coverage. No patient data are printed by checks.
stopifnot(!anyDuplicated(matched_pairs$patid),
          all(count(matched_pairs, ckd_patid)$n %in% 1:4))
checked <- matched_pairs %>%
  left_join(controls %>% select(patid, control_pracid = pracid,
                               regstartdate, gp_end_date), by = "patid") %>%
  left_join(cases %>% select(ckd_patid = patid, case_index = index_date,
                            case_pracid = pracid), by = "ckd_patid")

stopifnot(all(checked$pracid == checked$control_pracid),
          all(checked$pracid == checked$case_pracid),
          all(checked$index_date == checked$case_index),
          all(checked$regstartdate <= checked$index_date %m-% months(pre_index_months)),
          all(checked$gp_end_date >= checked$index_date %m+% months(post_index_months)))

stopifnot(
  all(
    rk_diabetes_advanced_ckd_matched$index_date >
      rk_diabetes_advanced_ckd_matched$dob %m+% years(18)
  ),
  all(
    rk_diabetes_egfr60_matched$index_date >
      rk_diabetes_egfr60_matched$dob %m+% years(18)
  )
)

# Verify dementia exclusions for BOTH groups at their shared case index date.
stopifnot(
  all(is.na(rk_diabetes_advanced_ckd_matched$first_dementia_date) |
        rk_diabetes_advanced_ckd_matched$first_dementia_date >
          rk_diabetes_advanced_ckd_matched$index_date %m+% months(dementia_lag_months)),
  all(is.na(rk_diabetes_egfr60_matched$first_dementia_date) |
        rk_diabetes_egfr60_matched$first_dementia_date >
          rk_diabetes_egfr60_matched$index_date %m+% months(dementia_lag_months))
)

# Inspect actual matching quality; nearest does not mean necessarily close.
pair_quality <- rk_diabetes_egfr60_matched %>%
  left_join(cases %>% select(ckd_patid = patid, case_dob = dob,
                            case_gender = gender,
                            case_regstartdate = regstartdate,
                            case_gp_end_date = gp_end_date), by = "ckd_patid") %>%
  summarise(
    median_age_gap_years = median(abs(as.numeric(dob - case_dob)) / 365.25),
    max_age_gap_years = if (n() > 0) max(abs(as.numeric(dob - case_dob)) / 365.25)
                       else NA_real_,
    proportion_same_gender = mean(gender == case_gender),
    median_regstart_gap_years =
      median(abs(as.numeric(regstartdate - case_regstartdate)) / 365.25),
    median_gpend_gap_years =
      median(abs(as.numeric(gp_end_date - case_gp_end_date)) / 365.25)
  )
print(pair_quality, width = Inf)

# Matching limitations:
# - Greedy matching need not maximise the number of retained cases.
# - Gender is approximate; diabetes type and duration are not matching variables.
# - max_age_gap_years remains Inf as in the uploaded script; set e.g. 5 if desired.
# - Review balance and account for differing n_controls in the outcome analysis.
# - This script does not construct dementia outcomes, censoring, or competing risks.

################################################################################
# Cache final matched cohorts AFTER eligibility and matching checks pass.
# These local outputs include the case-control links; a separate matched_pairs
# cache is unnecessary because the control table contains the same linkage.
rk_diabetes_advanced_ckd_matched_db <- cache_local_output(
  rk_diabetes_advanced_ckd_matched, "diabetes_advanced_ckd_matched",
  indexes = c("pracid", "index_date")
)
rk_diabetes_egfr60_matched_db <- cache_local_output(
  rk_diabetes_egfr60_matched, "diabetes_egfr60_matched",
  indexes = c("pracid", "ckd_patid", "index_date")
)

rk_diabetes_advanced_ckd_matched_db %>% count() %>% print()
rk_diabetes_egfr60_matched_db %>% count() %>% print()

# To load final results in a LATER session without rerunning this script:
# analysis <- cprd$analysis("ckd")
# rk_diabetes_advanced_ckd_matched_db <- rk_diabetes_advanced_ckd_matched_db %>%
#   analysis$cached("rk_adv_v1_diabetes_advanced_ckd_matched")
# rk_diabetes_egfr60_matched_db <- rk_diabetes_egfr60_matched_db %>%
#   analysis$cached("rk_adv_v1_diabetes_egfr60_matched")
#
# To resume before all final outputs exist, rerun the script with the same cache
# version. Completed SQL caches are reused; local R processing/matching reruns.
# If a server write was interrupted, verify that cache before trusting/reusing
# it. Existence alone cannot establish that an interrupted cache is complete.
