# Matched non-CKD controls for the diabetes CKD cohort
# Run AFTER the original cohort script in the same R session.
# No type 2 restriction is applied here. Rebuild the original CKD cache if it
# still contains only type 2 diabetes.
#
# Requirements: dplyr, tidyr, lubridate, MatchIt (current version).
# Matching is performed locally after collecting only the required columns.
# The outputs below are local R tables, NOT cached database tables.
#
# Assumptions requiring your review:
# - first_ckd_date is the intended CKD case index date. Change index_column if not.
# - The supplied dates are Date/POSIXt or ISO date strings, not numeric dates.
# - All patients must be registered 3 years pre-index date and at least 6 months post-index date.
# - Controls are absent from ckd_ids throughout the available records.
# - Membership in diabetes_cohort alone does not prove diabetes preceded index.
#   Apply diabetes-onset eligibility separately if required for the study.

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

index_column <- "first_ckd_date"

matching_columns <- c("patid", "pracid", "dob", "gender",
                      "regstartdate", "gp_enddate")

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
           across(any_of(c("dob", "regstartdate", "gp_enddate", "index_date")),
                  as.Date)) %>%
    drop_na(all_of(c("patid", "pracid", "dob", "gender",
                     "regstartdate", "gp_enddate"))) %>%
    filter(regstartdate <= gp_enddate) %>%
    arrange(patid)
}

cases <- prepare_patients(cases_raw) %>%
  drop_na(index_date) %>%
  mutate(window_start = index_date %m-% years(3),
         window_end = index_date %m+% years(0.5)) %>%
  filter(regstartdate <= window_start, gp_enddate >= window_end)

controls <- prepare_patients(controls_raw)
stopifnot(!any(cases$patid %in% controls$patid))

# Empty result with stable column types, including when no matches are found.
empty_pairs <- tibble(ckd_patid = character(), patid = character(),
                      index_date = as.Date(character()), pracid = character())

match_one_practice <- function(ca, co) {
  if (nrow(ca) == 0L || nrow(co) < 4L) return(empty_pairs)

  # Controls must span the PARTICULAR case's window. Rows are cases;
  # columns are controls. Different cases can have different index dates.
  eligible <- outer(as.numeric(ca$window_start),
                    as.numeric(co$regstartdate), `>=`) &
              outer(as.numeric(ca$window_end),
                    as.numeric(co$gp_enddate), `<=`)

  # Cases with fewer than four possible controls cannot form a complete set.
  keep_cases <- rowSums(eligible) >= 4L
  ca <- ca[keep_cases, , drop = FALSE]
  eligible <- eligible[keep_cases, , drop = FALSE]
  if (nrow(ca) == 0L) return(empty_pairs)

  keep_controls <- colSums(eligible) > 0L
  co <- co[keep_controls, , drop = FALSE]
  eligible <- eligible[, keep_controls, drop = FALSE]

  dat <- bind_rows(mutate(ca, CKD = 1L), mutate(co, CKD = 0L)) %>%
    mutate(gender = factor(gender),
           dob_number = as.numeric(dob),
           regstart_number = as.numeric(regstartdate),
           gpend_number = as.numeric(gp_enddate)) %>%
    as.data.frame()
  rownames(dat) <- dat$patid

  # These are approximate matching variables. Birth-date differences are
  # equivalent to age differences at the SAME case index date.
  # Drop practice-constant variables (including a single gender category).
  vars <- c("dob_number", "gender", "regstart_number", "gpend_number")
  vars <- vars[vapply(dat[vars], function(x) n_distinct(x) > 1L, logical(1))]
  form <- if (length(vars)) reformulate(vars, response = "CKD") else CKD ~ 1

  if (length(vars)) {
    distances <- MatchIt::mahalanobis_dist(form, data = dat)
  } else {
    distances <- matrix(0, nrow(ca), nrow(co))
  }
  dimnames(distances) <- list(ca$patid, co$patid)
  stopifnot(all(is.finite(distances)))

  # MatchIt treats Inf as a forbidden pairing.
  distances[!eligible] <- Inf

  fit <- matchit(
    form, data = dat,
    method = "nearest", distance = distances,
    exact = ~ pracid, ratio = 4, replace = FALSE, m.order = "random"
  )

  # ratio=4 requests four; it does not guarantee four are available.
  # Retain only complete sets, explicitly keeping the case-control linkage.
  mm <- fit$match.matrix
  mm <- mm[rowSums(!is.na(mm)) == 4L, , drop = FALSE]
  if (nrow(mm) == 0L) return(empty_pairs)

  tibble(ckd_patid = rep(rownames(mm), each = 4L),
         patid = as.vector(t(mm))) %>%
    left_join(ca %>% select(ckd_patid = patid, index_date, pracid),
              by = "ckd_patid")
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

matched_pairs <- bind_rows(empty_pairs, bind_rows(pair_list))

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
print(matching_summary)

# Structural checks: four unique controls, exact practice, shared case date,
# and control registration coverage. No patient data are printed by checks.
stopifnot(!anyDuplicated(matched_pairs$patid),
          all(count(matched_pairs, ckd_patid)$n == 4L))
checked <- matched_pairs %>%
  left_join(controls %>% select(patid, control_pracid = pracid,
                               regstartdate, gp_enddate), by = "patid") %>%
  left_join(cases %>% select(ckd_patid = patid, case_index = index_date,
                            case_pracid = pracid), by = "ckd_patid")
stopifnot(all(checked$pracid == checked$control_pracid),
          all(checked$pracid == checked$case_pracid),
          all(checked$index_date == checked$case_index),
          all(checked$regstartdate <= checked$index_date %m-% years(3)),
          all(checked$gp_enddate >= checked$index_date %m+% years(3)))

# Inspect actual matching quality; nearest does not mean necessarily close.
pair_quality <- rk_diabetes_nonckd_matched %>%
  left_join(cases %>% select(ckd_patid = patid, case_dob = dob,
                            case_gender = gender,
                            case_regstartdate = regstartdate,
                            case_gp_enddate = gp_enddate), by = "ckd_patid") %>%
  summarise(
    median_age_gap_years = median(abs(as.numeric(dob - case_dob)) / 365.25),
    max_age_gap_years = if (n() > 0) max(abs(as.numeric(dob - case_dob)) / 365.25)
                       else NA_real_,
    proportion_same_gender = mean(gender == case_gender),
    median_regstart_gap_years =
      median(abs(as.numeric(regstartdate - case_regstartdate)) / 365.25),
    median_gpend_gap_years =
      median(abs(as.numeric(gp_enddate - case_gp_enddate)) / 365.25)
  )
print(pair_quality)

# Limitations:
# - Greedy matching plus removing incomplete sets does NOT maximise the number
#   of complete 4:1 sets. Controls in incomplete sets are not reassigned here.
# - No maximum age/date gap is imposed; inspect pair_quality before analysis.
# - Gender and diabetes type are not exact-matching constraints.
# - Controls never meeting the CKD algorithm may still have unrecorded CKD.
# - Validate covariate balance in the retained cohort before outcome modelling.
