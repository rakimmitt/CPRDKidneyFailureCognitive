############################################################################################

#Setup
library(dplyr)
library(tidyr)
library(lubridate)
library(MatchIt)
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes_2024 = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix = "ckd"

###############################################################################################

# load ckd stages based on egfr only

analysis = cprd$analysis("all_patid")

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>% analysis$cached("ckd_stages_from_algorithm")

# load in acr data
clean_acr_medcodes <- clean_acr_medcodes %>%
  analysis$cached("clean_acr_medcodes", indexes=c("patid", "date", "testvalue"))

clean_acr_from_separate_medcodes <- clean_acr_from_separate_medcodes %>%
  analysis$cached("clean_acr_from_separate_medcodes", indexes=c("patid", "date", "testvalue"))

all_acr <- clean_acr_medcodes %>%
  select(patid, date, testvalue) %>%
  union_all(clean_acr_from_separate_medcodes %>%
              select(patid, date, testvalue))

# select those with acr >=3 mg/mmol
acr_high <- all_acr %>%
  filter(testvalue >= 3)

acr_span <- acr_high %>%
  group_by(patid) %>%
  summarise(
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    n_tests = n()
  )

# confirm 2 readings 3 months apart or longer
confirmed_acr3 <- acr_span %>%
  filter(n_tests >= 2 & datediff(max_date, min_date) >= 90) %>%
  mutate(confirmed_acr3_date = min_date) %>%
  select(patid, confirmed_acr3_date) %>%
  analysis$cached("confirmed_acr3", indexes = c("patid", "confirmed_acr3_date"))

# join with ckd stage
ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  left_join(confirmed_acr3, by = "patid") %>%
  mutate(
    stage_1 = case_when(
      is.na(stage_1) ~ sql("NULL"),
      is.na(confirmed_acr3_date) ~ sql("NULL"),
      confirmed_acr3_date <= stage_1 ~ stage_1,
      confirmed_acr3_date > stage_1  ~ confirmed_acr3_date
    ),
    stage_2 = case_when(
      is.na(stage_2) ~ sql("NULL"),
      is.na(confirmed_acr3_date) ~ sql("NULL"),
      confirmed_acr3_date <= stage_2 ~ stage_2,
      confirmed_acr3_date > stage_2  ~ confirmed_acr3_date
    )
  )  %>% 
  filter(!(is.na(stage_1) & is.na(stage_2) & is.na(stage_3a) & 
             is.na(stage_3b) & is.na(stage_4) & is.na(stage_5))) %>%
  analysis$cached("ckd_stages_from_algorithm_with_acr",
                  indexes = c("patid"))

# create list of ids of people with ckd and combine with type 2 diabetes cohort

analysis = cprd$analysis("diabetes_cohort")
practice_exclusion_ids <- practice_exclusion_ids %>% analysis$cached("practice_exclusion_ids")
gender_exclusion_ids <- gender_exclusion_ids %>% analysis$cached("gender_exclusion_ids")

ckd_ids <- ckd_stages_from_algorithm %>%
  mutate(
    first_ckd_date = as.Date(
      pmin(
        ifelse(is.na(stage_1), as.Date("2050-01-01"), stage_1),
        ifelse(is.na(stage_2), as.Date("2050-01-01"), stage_2),
        ifelse(is.na(stage_3a), as.Date("2050-01-01"), stage_3a),
        ifelse(is.na(stage_3b), as.Date("2050-01-01"), stage_3b),
        ifelse(is.na(stage_4), as.Date("2050-01-01"), stage_4),
        ifelse(is.na(stage_5), as.Date("2050-01-01"), stage_5),
        na.rm = TRUE
      )
    )
  ) %>%
  mutate(first_ckd_date = ifelse(first_ckd_date == as.Date("2050-01-01"), NA, first_ckd_date)) %>%
  group_by(patid) %>%
  dbplyr::window_order(first_ckd_date) %>%
  distinct(patid, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-contains("stage"), -confirmed_acr3_date) %>% 
  anti_join(practice_exclusion_ids, by="patid") %>% anti_join(gender_exclusion_ids, by="patid") %>% 
  analysis$cached("ckd_ids", unique_indexes="patid")

ckd_ids %>% count() # 913743

# load diabetes cohort
analysis = cprd$analysis("all")

diabetes_cohort <- diabetes_cohort %>% 
  analysis$cached("diabetes_cohort")

rk_diabetes_ckd_cohort <- diabetes_cohort %>%
  inner_join(ckd_ids, by = "patid") %>%
  analysis$cached("rk_diabetes_ckd_cohort", unique_indexes = "patid")

rk_diabetes_ckd_cohort  %>% count()

# Matched non-CKD controls for the diabetes CKD cohort
# Matching is performed locally after collecting only the required columns.
# The outputs below are local R tables, NOT cached database tables.
#
# Assumptions:
# - first_ckd_date is the intended CKD case index date. Change index_column if not.
# - All patients must be registered 3 years pre-index date and at least 6 months post-index date.
# - Controls are absent from ckd_ids throughout the available records.
# - Membership in diabetes_cohort alone does not prove diabetes preceded index - need to apply diabetes-onset eligibility separately if required

index_column <- "first_ckd_date"

matching_columns <- c("patid", "pracid", "dob", "gender",
                      "regstartdate", "gp_end_date")

cases_raw <- rk_diabetes_ckd_cohort %>%
  select(all_of(c(matching_columns, index_column))) %>%
  collect() %>%
  rename(index_date = all_of(index_column))

controls_raw <- diabetes_cohort %>%
  anti_join(ckd_ids, by = "patid") %>%
  anti_join(practice_exclusion_ids, by = "patid") %>%
  anti_join(gender_exclusion_ids, by = "patid") %>%
  select(all_of(matching_columns)) %>%
  collect()

# Do not silently discard repeated patients or conflicting demographics.
stopifnot(!anyDuplicated(cases_raw$patid),
          !anyDuplicated(controls_raw$patid))

prepare_patients <- function(x) {
  x %>%
    mutate(patid = as.character(patid), pracid = as.character(pracid),
           gender = as.character(gender),
           across(any_of(c("dob", "regstartdate", "gp_end_date", "index_date")),
                  as.Date)) %>%
    drop_na(all_of(c("patid", "pracid", "dob", "gender",
                     "regstartdate", "gp_end_date"))) %>%
    filter(regstartdate <= gp_end_date) %>%
    arrange(patid)
}

cases <- prepare_patients(cases_raw) %>%
  drop_na(index_date) %>%
  mutate(
    window_start = index_date %m-% months(6),
    window_end = index_date %m+% months(6)
  ) %>%
  filter(
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

rk_diabetes_nonckd_matched <- matched_pairs %>%
  left_join(controls %>% select(-pracid), by = "patid")

rk_diabetes_ckd_matched <- cases %>%
  semi_join(matched_pairs, by = c("patid" = "ckd_patid"))

matching_summary <- tibble(
  original_ckd_patients = nrow(cases_raw),
  ckd_patients_with_complete_data_and_window = nrow(cases),
  matched_ckd_patients = nrow(rk_diabetes_ckd_matched),
  matched_controls = nrow(rk_diabetes_nonckd_matched)
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
          all(checked$regstartdate <= checked$index_date %m-% months(6)),
          all(checked$gp_end_date >= checked$index_date %m+% months(6)))

stopifnot(
  all(
    rk_diabetes_ckd_matched$index_date >
      rk_diabetes_ckd_matched$dob %m+% years(18)
  ),
  all(
    rk_diabetes_nonckd_matched$index_date >
      rk_diabetes_nonckd_matched$dob %m+% years(18)
  )
)

# Inspect actual matching quality; nearest does not mean necessarily close.
pair_quality <- rk_diabetes_nonckd_matched %>%
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

# Limitations:
# - Greedy matching plus removing incomplete sets does NOT maximise the number
#   of complete 4:1 sets. Controls in incomplete sets are not reassigned here.
# - No maximum age/date gap is imposed; inspect pair_quality before analysis.
# - Controls never meeting the CKD algorithm may still have unrecorded CKD.
# - Validate covariate balance in the retained cohort before outcome modelling.