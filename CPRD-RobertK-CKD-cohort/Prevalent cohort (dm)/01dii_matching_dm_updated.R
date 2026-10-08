# Match advanced CKD cases to date-eligible controls, with replacement.

# Input: the broad matching_pool from 01c (one row per patient).
# Output: one case row plus up to four control rows per matched set.

# A patient may occur in multiple sets and may have different roles over time.

# Setup

library(tidyverse)
library(lubridate)
library(aurum)
library(MatchIt)

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")

analysis = cprd$analysis("rk_ckd")

matching_pool <- matching_pool %>%
  analysis$cached("matching_pool")

age_penalty <- 0 # 0 reproduces original distance
match_with_replacement <- TRUE
max_controls <- 4L
max_age_gap_years <- Inf
min_registration_months <- 6L
dementia_exclusion_months <- 3L
case_date_cutoff <- as.Date("2008-01-01")
set.seed(123)

stopifnot(max_controls >= 1L, max_controls == as.integer(max_controls),
          max_age_gap_years >= 0)

columns <- c(
  "patid", "pracid", "dob", "gender", "ethnicity_5cat", "imd_decile",
  "regstartdate", "gp_end_date", "hes_end_date", "death_date", "with_hes",
  "dm_diag_date_all", "earliest_all_dementia",
  "ckd_stage_3_start_date", "advanced_ckd_index_date"
)
missing_columns <- setdiff(columns, colnames(matching_pool))
if (length(missing_columns)) {
  stop("Missing columns in matching_pool: ", paste(missing_columns, collapse = ", "))
}

# These two exclusions do not depend on a case date.
pool_raw <- matching_pool %>%
  filter(with_hes == 1L, !is.na(dm_diag_date_all)) %>%
  select(all_of(columns)) %>%
  mutate(patid = as.character(patid), pracid = as.character(pracid)) %>%
  collect()
stopifnot(!anyNA(pool_raw$patid), !anyDuplicated(pool_raw$patid))

# 2. Prepare the broad pool and select eligible cases
event_day <- function(x) coalesce(as.numeric(x), Inf)
date_columns <- c(
  "dob", "regstartdate", "gp_end_date", "hes_end_date", "death_date",
  "dm_diag_date_all", "earliest_all_dementia",
  "ckd_stage_3_start_date", "advanced_ckd_index_date"
)

pool <- pool_raw %>%
  mutate(
    across(all_of(date_columns), as.Date),
    gender = as.character(gender),
    ethnicity_5cat = coalesce(as.character(ethnicity_5cat), "Missing"),
    imd_decile = as.numeric(na_if(as.character(imd_decile), "Missing")),
    imd_missing = as.integer(is.na(imd_decile)),
    imd_for_matching = coalesce(imd_decile, 5.5),
    adult_date = dob %m+% years(18),
    dementia_day = event_day(earliest_all_dementia),
    diabetes_day = as.numeric(dm_diag_date_all),
    control_stop_day = pmin(event_day(ckd_stage_3_start_date),
                            event_day(advanced_ckd_index_date)),
    # GP and HES end dates must both be known. An absent death date does
    # not shorten observation; a recorded death does.
    followup_end_day = pmin(as.numeric(gp_end_date),
                           as.numeric(hes_end_date), event_day(death_date))
  ) %>%
  filter(
    !is.na(dm_diag_date_all),
    gp_end_date >= regstartdate %m+% months(min_registration_months)
  ) %>%
  arrange(patid)

# Identify cases (advanced_ckd_index_date after 1/1/2008, adult, no dementia within 3 months of index)
cases <- pool %>%
  filter(!is.na(advanced_ckd_index_date)) %>%
  mutate(index_date = advanced_ckd_index_date) %>%
  filter(
    index_date >= case_date_cutoff,
    adult_date <= index_date,
    regstartdate <= index_date,
    as.numeric(index_date) <= followup_end_day,
    dementia_day > as.numeric(index_date %m+% months(dementia_exclusion_months))
  )

# Do not exclude future cases or future dementia/CKD diagnoses globally.
# Their eligibility will be recalculated at each prospective case's date.

# Controls must have eligible information in order to be matched
controls <- pool %>%
  drop_na(pracid, dob, gender, regstartdate, followup_end_day)

stopifnot(!anyDuplicated(cases$patid), !anyDuplicated(controls$patid))

empty_pairs <- tibble(matched_case_patid = character(), patid = character())

# Match within a practice
match_practice <- function(ca, co) {
  ca <- ca %>% drop_na(pracid, dob, gender, index_date)
  if (!nrow(ca) || is.null(co) || !nrow(co)) return(empty_pairs)

  # outer(x, y, comparison) tests EVERY element of x against EVERY element
  # of y. Rows below are cases, columns are potential controls.
  # One control can therefore be eligible in one row and ineligible in another.
  
  # define dementia and diabetes status at case index date
  case_day <- as.numeric(ca$index_date)
  dementia_cutoff <- as.numeric(
    ca$index_date %m+% months(dementia_exclusion_months)
  )
  case_dm <- ca$diabetes_day <= case_day
  control_dm <- outer(case_day, co$diabetes_day, `>=`)
  same_dm <- sweep(control_dm, 1, case_dm, `==`)

  eligible <-
    outer(ca$patid, co$patid, `!=`) &                    # no self-match
    outer(case_day, as.numeric(co$regstartdate), `>=`) & # already registered
    outer(case_day, co$followup_end_day, `<=`) &         # observable at index
    outer(case_day, as.numeric(co$adult_date), `>=`) &   # at least 18
    outer(dementia_cutoff, co$dementia_day, `<`) &       # no early dementia
    outer(case_day, co$control_stop_day, `<`) &          # before renal cut-off
    same_dm &                                            # exact diabetes status
    (abs(outer(as.numeric(ca$dob), as.numeric(co$dob), `-`)) <=
       max_age_gap_years * 365.25)
  eligible[is.na(eligible)] <- FALSE

  # Remove impossible rows/columns ONLY from this local matching call.
  keep_cases <- rowSums(eligible) > 0L
  ca <- ca[keep_cases, , drop = FALSE]
  eligible <- eligible[keep_cases, , drop = FALSE]
  if (!nrow(ca)) return(empty_pairs)
  keep_controls <- colSums(eligible) > 0L
  co <- co[keep_controls, , drop = FALSE]
  eligible <- eligible[, keep_controls, drop = FALSE]

  # Unique role IDs are required because the same patid may be in both pools.
  case_keys <- paste0("case_", ca$patid)
  control_keys <- paste0("control_", co$patid)
  dat <- bind_rows(mutate(ca, is_case = 1L), mutate(co, is_case = 0L)) %>%
    mutate(dob_days = as.numeric(dob),
           regstart_days = as.numeric(regstartdate),
           gender = factor(gender), ethnicity_5cat = factor(ethnicity_5cat)) %>%
    as.data.frame()
  rownames(dat) <- c(case_keys, control_keys)

  # Soft matching: nearer overall Mahalanobis distance is preferred.
  # Missing IMD remains NA in outputs; 5.5 is only a matching placeholder.
  # Its separate indicator allows missingness to contribute to the distance.
  variables <- c("dob_days", "gender", "ethnicity_5cat",
                 "imd_for_matching", "imd_missing", "regstart_days")
  stopifnot(all(complete.cases(dat[variables])))
  variables <- variables[
    vapply(dat[variables], function(x) n_distinct(x) > 1L, logical(1))
  ]
  form <- if (length(variables)) reformulate(variables, response = "is_case") else is_case ~ 1
  distances <- if (nrow(ca) == 1L && nrow(co) == 1L || !length(variables)) {
    matrix(0, nrow(ca), nrow(co))
  } else {
    MatchIt::mahalanobis_dist(form, data = dat)
  }

# Calculate age differences for the remaining cases and controls.
age_gap_years <- abs(
  outer(as.numeric(ca$dob), as.numeric(co$dob), `-`)
) / 365.25

# Retain the existing matching distance, but penalise age gaps further.
distances <- distances + age_penalty * age_gap_years

  dimnames(distances) <- list(case_keys, control_keys)
  stopifnot(all(is.finite(distances)))
  distances[!eligible] <- Inf

  fit <- MatchIt::matchit(
    form, data = dat, method = "nearest", distance = distances,
    ratio = min(max_controls, nrow(co)), replace = match_with_replacement, m.order = "random"
  )
  mm <- fit$match.matrix
  if (is.null(mm) || !length(mm)) return(empty_pairs)
  tibble(case_key = rep(rownames(mm), each = ncol(mm)),
         control_key = as.vector(t(mm))) %>%
    filter(!is.na(control_key)) %>%
    transmute(matched_case_patid = ca$patid[match(case_key, case_keys)],
              patid = co$patid[match(control_key, control_keys)])
}

# Each practice has its own distance matrix, enforcing exact pracid.
# Very large practices can still require considerable RAM.
case_groups <- split(seq_len(nrow(cases)), cases$pracid)
control_groups <- split(seq_len(nrow(controls)), controls$pracid)
pair_list <- vector("list", length(case_groups))
for (i in seq_along(case_groups)) {
  practice <- names(case_groups)[i]
  ca <- cases[case_groups[[i]], , drop = FALSE]
  rows <- control_groups[[practice]]
  co <- controls[if (is.null(rows)) integer() else rows, , drop = FALSE]
  pair_list[[i]] <- match_practice(ca, co)
  if (i %% 100L == 0L) message("Processed ", i, " practices")
}
pairs <- bind_rows(empty_pairs, bind_rows(pair_list))

#   Outputs: preserve every matched episode 
set_sizes <- pairs %>% count(matched_case_patid, name = "n_controls")

# Keep only cases with at least one matched control.
matched_cases <- cases %>%
  inner_join(set_sizes, by = c("patid" = "matched_case_patid"))

matched_controls <- pairs %>%
  left_join(controls, by = "patid") %>%
  left_join(cases %>% transmute(matched_case_patid = patid,
                               matched_case_index_date = index_date),
            by = "matched_case_patid") %>%
  left_join(set_sizes, by = "matched_case_patid") %>%
  mutate(index_date = matched_case_index_date)

matched_cohort <- bind_rows(
  matched_cases %>% mutate(is_case = 1L, matched_case_patid = patid,
                          matched_case_index_date = index_date),
  matched_controls %>% mutate(is_case = 0L)
) %>%
  mutate(
    match_record_id = row_number(),
    dm_at_index = as.integer(dm_diag_date_all <= index_date),
    # Descriptive set weights: controls in each set sum to one.
    weight = if_else(is_case == 1L, 1, 1 / n_controls)
  ) %>%
  select(match_record_id, patid, is_case, matched_case_patid,
         matched_case_index_date, index_date, n_controls, weight, dm_at_index,
         all_of(setdiff(columns, "patid")))

#   Independent pair checks BEFORE caching 
stopifnot(
  !anyNA(pairs),
  !anyDuplicated(pairs[c("matched_case_patid", "patid")]),
  all(pairs$patid != pairs$matched_case_patid),
  all(matched_cases$n_controls <= max_controls),
  nrow(matched_cases) == n_distinct(pairs$matched_case_patid),
  all(matched_cases$n_controls >= 1L),
  nrow(matched_controls) == sum(matched_cases$n_controls)
)

if (!match_with_replacement) {
  stopifnot(!anyDuplicated(matched_controls$patid))
}

pair_comparison <- matched_controls %>%
  left_join(cases %>% transmute(
    matched_case_patid = patid, case_index = index_date, case_pracid = pracid,
    case_dob = dob, case_gender = gender, case_ethnicity = ethnicity_5cat,
    case_imd = imd_decile, case_regstart = regstartdate,
    case_dm_date = dm_diag_date_all
  ), by = "matched_case_patid") %>%
  mutate(
    age_gap_years = abs(as.numeric(dob - case_dob)) / 365.25,
    same_gender = gender == case_gender,
    same_ethnicity = ethnicity_5cat == case_ethnicity,
    imd_gap = abs(imd_decile - case_imd),
    both_imd_observed = !is.na(imd_decile) & !is.na(case_imd),
    same_imd_missingness = is.na(imd_decile) == is.na(case_imd),
    registration_gap_years = abs(as.numeric(regstartdate - case_regstart)) / 365.25,
    same_diabetes_status = (dm_diag_date_all <= case_index) ==
      (case_dm_date <= case_index)
  )
stopifnot(
  all(pair_comparison$pracid == pair_comparison$case_pracid),
  all(pair_comparison$index_date == pair_comparison$case_index),
  all(pair_comparison$matched_case_index_date == pair_comparison$case_index),
  all(pair_comparison$adult_date <= pair_comparison$case_index),
  all(pair_comparison$regstartdate <= pair_comparison$case_index),
  all(as.numeric(pair_comparison$case_index) <= pair_comparison$followup_end_day),
  all(pair_comparison$gp_end_date >=
        pair_comparison$regstartdate %m+% months(min_registration_months)),
  all(pair_comparison$dementia_day > as.numeric(
    pair_comparison$case_index %m+% months(dementia_exclusion_months))),
  all(as.numeric(pair_comparison$case_index) < pair_comparison$control_stop_day),
  all(pair_comparison$same_diabetes_status),
  all(pair_comparison$age_gap_years <= max_age_gap_years)
)

matching_summary <- tibble(
  pool_after_hes_and_known_diabetes = nrow(pool_raw),
  pool_after_registration_rule = nrow(pool),
  eligible_cases = nrow(cases),
  cases_with_controls = nrow(matched_cases),
  cases_without_controls = nrow(cases) - nrow(matched_cases),
  control_records = nrow(matched_controls),
  unique_control_patients = n_distinct(matched_controls$patid),
  patients_serving_as_both_case_and_control =
    length(intersect(matched_cases$patid, matched_controls$patid))
)
print(matching_summary, width = Inf)
matched_cases %>% count(n_controls, name = "number_of_cases") %>% print(n = Inf)
control_reuse <- matched_controls %>% count(patid, name = "times_used")
control_reuse %>% count(times_used, name = "number_of_patients") %>% print(n = Inf)

# Save the complete episode table
# patid is NOT unique now. Keep the local table and the database reference
# in separate objects. This replaces the named output table when run.

if (nrow(matched_cohort) == 0L) {
  stop("No matched sets were produced; the existing output table was not replaced.")
}

matched_cohort_db <- copy_to(
  dest = analysis$.con, df = matched_cohort,
  name = dbplyr::in_schema(analysis$.analysisDb, "rk_ckd_matched_cohort"),
  overwrite = TRUE, temporary = FALSE,
  unique_indexes = "match_record_id",
  indexes = c("patid", "is_case", "matched_case_patid", "index_date")
)
# Wait for copy_to() to finish; an interrupted client upload may be incomplete.

# Balance summary statistics

# Balance summaries that preserve reused control episodes ----------------
# Gross means unweighted AFTER matching, not balance in the original pool.
# Do not deduplicate control patid: age/diabetes status can vary by episode.
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_wmean <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) NA_real_ else weighted.mean(x[ok], w[ok])
}
safe_wmedian <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  # Lower weighted median: first value reaching half of total weight.
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

summarise_balance <- function(include_unmatched_cases = FALSE) {
  b <- matched_cohort %>%
    mutate(group = if_else(is_case == 1L, "Cases", "Controls"),
           age = as.numeric(index_date - dob) / 365.25,
           # Duration before index compares registration starts at shared dates.
           registration_years_at_index = as.numeric(index_date - regstartdate) / 365.25,
           gender = coalesce(as.character(gender), "Missing"))
  if (!nrow(b) || !all(c("Cases", "Controls") %in% b$group)) {
    message("Both groups are required for balance summaries.")
    return(NULL)
  }
  numeric_summary <- b %>%
    pivot_longer(c(age, imd_decile, registration_years_at_index),
                 names_to = "variable", values_to = "value") %>%
    group_by(variable, group) %>%
    summarise(n = n(), missing_n = sum(is.na(value)),
              missing_percent = 100 * mean(is.na(value)),
              weighted_missing_percent = 100 * weighted.mean(is.na(value), weight),
              mean = safe_mean(value), sd = sd(value, na.rm = TRUE),
              median = median(value, na.rm = TRUE),
              weighted_mean = safe_wmean(value, weight),
              weighted_median = safe_wmedian(value, weight), .groups = "drop")
  categorical_summary <- b %>%
    mutate(imd_category = coalesce(as.character(imd_decile), "Missing"),
           imd_missing = if_else(is.na(imd_decile), "Missing", "Observed"),
           diabetes_status = if_else(dm_at_index == 1L, "Diagnosed", "Not yet diagnosed")) %>%
    pivot_longer(c(gender, ethnicity_5cat, imd_category, imd_missing, diabetes_status),
                 names_to = "variable", values_to = "category") %>%
    group_by(variable, category, group) %>%
    summarise(n = n(), weighted_n = sum(weight), .groups = "drop") %>%
    complete(nesting(variable, category), group = c("Cases", "Controls"),
             fill = list(n = 0L, weighted_n = 0)) %>%
    group_by(variable, group) %>%
    mutate(percent = 100 * n / sum(n),
           weighted_percent = 100 * weighted_n / sum(weighted_n)) %>% ungroup()

  # Fixed case SD for numeric variables; binary SD for each category.
  # IMD numeric balance is conditional on observation; also assess missingness.
  numeric_smd <- numeric_summary %>%
    select(variable, group, mean, weighted_mean, sd) %>%
    pivot_wider(names_from = group, values_from = c(mean, weighted_mean, sd)) %>%
    transmute(variable, category = NA_character_, denominator = sd_Cases,
              unweighted_difference = mean_Cases - mean_Controls,
              weighted_difference = weighted_mean_Cases - weighted_mean_Controls)
  categorical_smd <- categorical_summary %>%
    select(variable, category, group, percent, weighted_percent) %>%
    pivot_wider(names_from = group, values_from = c(percent, weighted_percent)) %>%
    transmute(variable, category,
              denominator = sqrt(percent_Cases / 100 * (1 - percent_Cases / 100)),
              unweighted_difference = (percent_Cases - percent_Controls) / 100,
              weighted_difference = (weighted_percent_Cases - weighted_percent_Controls) / 100)
  smd_summary <- bind_rows(numeric_smd, categorical_smd) %>%
    mutate(unweighted_smd = if_else(denominator > 0, unweighted_difference / denominator, NA_real_),
           weighted_smd = if_else(denominator > 0, weighted_difference / denominator, NA_real_),
           absolute_weighted_smd = abs(weighted_smd))
  print(numeric_summary, n = Inf, width = Inf)
  print(categorical_summary, n = Inf, width = Inf)
  print(smd_summary, n = Inf, width = Inf)
  list(numeric_summary = numeric_summary, categorical_summary = categorical_summary,
       smd_summary = smd_summary)
}

message("Balance: cases with at least one control")
balance_matched <- summarise_balance(FALSE)

# Pair summaries: each pair counts equally; larger sets contribute more.
pair_quality <- pair_comparison %>% summarise(
  n_pairs = n(),
  mean_age_gap = safe_mean(age_gap_years),
  median_age_gap = median(age_gap_years, na.rm = TRUE),
  p95_age_gap = as.numeric(quantile(age_gap_years, 0.95, na.rm = TRUE)),
  largest_age_gap = safe_max(age_gap_years),
  percent_same_gender = 100 * safe_mean(same_gender),
  percent_same_ethnicity = 100 * safe_mean(same_ethnicity),
  percent_same_diabetes_status = 100 * safe_mean(same_diabetes_status),
  percent_same_imd_missingness = 100 * safe_mean(same_imd_missingness),
  n_pairs_with_both_imd_observed = sum(both_imd_observed),
  mean_imd_gap_observed = safe_mean(imd_gap),
  largest_imd_gap_observed = safe_max(imd_gap),
  mean_registration_gap_years = safe_mean(registration_gap_years)
)
print(pair_quality, width = Inf)

# Notes for subsequent analysis:
# - Join patient-level tables to these records without distinct(patid).
# - Derive time-dependent characteristics using each row's assigned index_date.
# - is_case describes this episode, not a permanent patient characteristic.
# - Repeated patients are not independent observations for outcome modelling.
# - Control eligibility is checked at index. This script does not implement
#   subsequent censoring at CKD onset or an outcome-analysis time origin.
# - Excluding recorded dementia through index + 3 months does not guarantee
#   three months of observable follow-up. No such follow-up minimum is imposed.
# - SMDs are descriptive. Zero case variance gives NA; inspect raw differences.
# - These weights describe matched sets and are not a complete survival model.
