#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(tidyverse)
library(MatchIt)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("rk_ckd")

advanced_ckd_cohort <- advanced_ckd_cohort %>%
  analysis$cached("advanced_ckd_cohort")
  advanced_ckd_cohort %>% count()

  non_ckd_cohort <- non_ckd_cohort %>%
  analysis$cached("non_ckd_cohort")
  non_ckd_cohort %>% count()

# Settings
max_controls <- 4L
max_age_gap_years <- 3
set.seed(123)

# 1. Load the required columns; both groups must have HES linkage
columns <- c(
  "patid", "pracid", "dob", "gender", "regstartdate",
  "gp_end_date", "hes_end_date", "with_hes", "ethnicity_5cat", "imd_decile"
)

cases <- advanced_ckd_cohort %>%
  filter(with_hes == 1) %>%
  select(all_of(c(columns, "index_date"))) %>%
  collect()

controls <- non_ckd_cohort %>%
  filter(with_hes == 1) %>%
  select(all_of(columns)) %>%
  collect()

# Preserve patient IDs as text and ensure dates have Date class
prepare <- function(x) {
  x %>%
    mutate(
      across(c(patid, pracid, gender), as.character),

      # Treat missing ethnicity as an explicit matching category.
      ethnicity_5cat = coalesce(as.character(ethnicity_5cat), "Missing"),

      # Preserve the original IMD, including NA.
      # na_if() also handles "Missing" if previously assigned.
      imd_decile = as.numeric(
      na_if(as.character(imd_decile), "Missing")
      ),

      # Additional variables used only for matching.
      imd_missing = as.integer(is.na(imd_decile)),
      imd_for_matching = coalesce(imd_decile, 5.5),

      across(
        any_of(c("dob", "regstartdate", "gp_end_date",
                 "hes_end_date", "index_date")),
        as.Date
      )
    ) %>%
    arrange(patid)
}

cases <- prepare(cases)
controls <- prepare(controls)

# Each patient must appear once, and cases cannot also be controls
stopifnot(
  !anyNA(cases$patid), !anyNA(controls$patid),
  !anyDuplicated(cases$patid), !anyDuplicated(controls$patid),
  !any(cases$patid %in% controls$patid)
)

# Controls with missing matching/eligibility information cannot be matched
controls <- controls %>%
  drop_na(patid, pracid, dob, gender,
          regstartdate, gp_end_date, hes_end_date)

empty_pairs <- tibble(
  matched_case_patid = character(),
  patid = character()
)

# 2. Match within one practice
match_practice <- function(ca, co) {

  # Cases with missing information remain in the final case table,
  # but cannot receive controls.
  ca <- ca %>% drop_na(pracid, dob, gender, index_date)

  if (nrow(ca) == 0L || is.null(co) || nrow(co) == 0L) {
    return(empty_pairs)
  }

  # Rows = cases; columns = controls.
  # All dates are inclusive. Age gap uses 365.25 days per year.
  eligible <-
    outer(as.numeric(ca$index_date),
          as.numeric(co$regstartdate), `>=`) &
    outer(as.numeric(ca$index_date),
          as.numeric(co$gp_end_date), `<=`) &
    outer(as.numeric(ca$index_date),
          as.numeric(co$hes_end_date), `<=`) &
    abs(outer(as.numeric(ca$dob), as.numeric(co$dob), `-`)) <=
      max_age_gap_years * 365.25

  # Remove patients with no possible pairing from this matching call only
  keep_cases <- rowSums(eligible) > 0L
  ca <- ca[keep_cases, , drop = FALSE]
  eligible <- eligible[keep_cases, , drop = FALSE]

  if (nrow(ca) == 0L) return(empty_pairs)

  keep_controls <- colSums(eligible) > 0L
  co <- co[keep_controls, , drop = FALSE]
  eligible <- eligible[, keep_controls, drop = FALSE]

  dat <- bind_rows(
    mutate(ca, is_case = 1L),
    mutate(co, is_case = 0L)
  ) %>%
    mutate(
      dob_days = as.numeric(dob),
      gender = factor(gender),
      ethnicity_5cat = factor(ethnicity_5cat)
    ) %>%
    as.data.frame()

  rownames(dat) <- dat$patid

  # Omit variables that are constant within this practice
  variables <- c("dob_days", "gender", "ethnicity_5cat", "imd_for_matching", "imd_missing")
  variables <- variables[
    vapply(dat[variables], function(x) n_distinct(x) > 1L, logical(1))
  ]

  formula <- if (length(variables)) {
    reformulate(variables, response = "is_case")
  } else {
    is_case ~ 1
  }

  # A single possible pair needs no distance ranking
  distances <- if (nrow(ca) == 1L && nrow(co) == 1L) {
    matrix(0, 1L, 1L)
  } else if (length(variables)) {
    MatchIt::mahalanobis_dist(formula, data = dat)
  } else {
    matrix(0, nrow(ca), nrow(co))
  }

  dimnames(distances) <- list(ca$patid, co$patid)
  stopifnot(all(is.finite(distances)))

  # Prohibit pairs that fail the date or age-gap requirements
  distances[!eligible] <- Inf

  fit <- matchit(
    formula,
    data = dat,
    method = "nearest",
    distance = distances,
    ratio = min(max_controls, nrow(co)),
    replace = FALSE,
    m.order = "random"
  )

  # Extract actual matches; NA entries are unfilled control slots
  mm <- fit$match.matrix

  tibble(
    matched_case_patid = rep(rownames(mm), each = ncol(mm)),
    patid = as.vector(t(mm))
  ) %>%
    filter(!is.na(patid))
}

# 3. Run separately within each practice: this guarantees exact pracid
case_groups <- split(cases, cases$pracid)
control_groups <- split(controls, controls$pracid)

pairs <- lapply(names(case_groups), function(practice) {
  match_practice(case_groups[[practice]], control_groups[[practice]])
}) %>%
  bind_rows(empty_pairs)

# 4. Create the matched-control table with the requested case linkage
matched_controls <- pairs %>%
  left_join(controls, by = "patid") %>%
  left_join(
    cases %>%
      transmute(
        matched_case_patid = patid,
        matched_case_index_date = index_date,
        index_date = index_date
      ),
    by = "matched_case_patid"
  )

# Retain ALL HES-linked cases, including cases with zero controls
matched_cases <- cases %>%
  left_join(
    pairs %>% count(matched_case_patid, name = "n_controls"),
    by = c("patid" = "matched_case_patid")
  ) %>%
  mutate(n_controls = coalesce(n_controls, 0L))

# Optional combined table: case rows carry their own ID and index date
matched_cohort <- bind_rows(
  matched_cases %>%
    mutate(
      is_case = 1L,
      matched_case_patid = patid,
      matched_case_index_date = index_date
    ),
  matched_controls %>% mutate(is_case = 0L)
)

analysis = cprd$analysis("rk_ckd")

matched_cohort <- copy_to(
  dest = analysis$.con,
  df = matched_cohort %>%
    select(patid, is_case, matched_case_patid, matched_case_index_date,
           index_date, dob, gender, ethnicity_5cat, imd_decile, regstartdate,
           gp_end_date, hes_end_date),
  name = dbplyr::in_schema(
  analysis$.analysisDb,
  "rk_ckd_matched_cohort"
  ),
  overwrite = TRUE,
  temporary = FALSE,
  unique_indexes = "patid",
  indexes = c("is_case", "matched_case_patid", "index_date")
)


# 5. Check control reuse and display the number of controls per case
stopifnot(
  !anyDuplicated(matched_controls$patid),
  all(matched_cases$n_controls <= max_controls),
  nrow(matched_cases) == nrow(cases)
)

matched_cases %>%
  count(n_controls, name = "number_of_cases") %>%
  print(n = Inf)






########

# Checking quality of matches

# Combine matched cases and controls
# NB because of the variable ratio matching, there may be differences in overall age and gender distributions
# Pair-wise matching should be tighter, and we can summarise weighted matching characteristics

# Balance and pair-quality summaries only
# Inputs: local matched_cases and matched_controls tables.
# No matching, database writes or changes to these input tables are performed.
# Gross = unweighted balance AFTER matching, not pre-match balance.

library(tidyverse)

# 1. Prepare separate copies for summaries and count actual controls per case.
# This avoids relying on stale n_controls or any previous weight columns.
balance_cases <- matched_cases %>%
  select(patid, dob, gender, ethnicity_5cat, imd_decile, index_date) %>%
  mutate(patid = as.character(patid))

balance_controls <- matched_controls %>%
  select(patid, matched_case_patid, dob, gender, ethnicity_5cat,
         imd_decile, matched_case_index_date) %>%
  mutate(across(c(patid, matched_case_patid), as.character))

stopifnot(
  !anyNA(balance_cases$patid), !anyNA(balance_controls$patid),
  !anyDuplicated(balance_cases$patid),
  !anyDuplicated(balance_controls$patid),
  !any(balance_cases$patid %in% balance_controls$patid),
  all(balance_controls$matched_case_patid %in% balance_cases$patid)
)

set_sizes <- balance_controls %>%
  count(matched_case_patid, name = "n_controls")

balance_cases <- balance_cases %>%
  left_join(set_sizes, by = c("patid" = "matched_case_patid")) %>%
  mutate(n_controls = coalesce(n_controls, 0L))

matching_counts <- balance_cases %>%
  count(n_controls, name = "number_of_cases")
print(matching_counts, n = Inf)

# Include ALL cases, including those with zero controls.
balance_data <- bind_rows(
  balance_cases %>%
    transmute(
      group = "Cases",
      dob, gender, ethnicity_5cat, imd_decile,
      index_date,
      weight = 1
    ),

  balance_controls %>%
    left_join(set_sizes, by = "matched_case_patid") %>%
    transmute(
      group = "Controls",
      dob, gender, ethnicity_5cat, imd_decile,
      index_date = matched_case_index_date,
      weight = 1 / n_controls
    )
) %>%
  mutate(
    age = as.numeric(as.Date(index_date) - as.Date(dob)) / 365.25,
    gender = coalesce(as.character(gender), "Missing"),
    ethnicity_5cat = coalesce(as.character(ethnicity_5cat), "Missing"),
    imd_decile = as.numeric(na_if(as.character(imd_decile), "Missing"))
  ) %>%
  select(group, age, gender, ethnicity_5cat, imd_decile, weight)


stopifnot(
  all(is.finite(balance_data$weight)),
  all(balance_data$weight > 0)
)

# Helpers return NA, rather than Inf/NaN, if nothing is observed.
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_weighted_mean <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  if (!any(keep)) NA_real_ else weighted.mean(x[keep], w[keep])
}

# 2. Population summaries and SMDs.
if (nrow(balance_data) == 0L) {
  message("No matches found: population balance summaries are unavailable.")
} else {
  numeric_summary <- balance_data %>%
    pivot_longer(c(age, imd_decile), names_to = "variable", values_to = "value") %>%
    group_by(variable, group) %>%
    summarise(
      n = n(), observed_n = sum(!is.na(value)),
      missing_n = sum(is.na(value)),
      missing_percent = 100 * mean(is.na(value)),
      weighted_missing_percent = 100 * weighted.mean(is.na(value), weight),
      mean = safe_mean(value), sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      weighted_mean = safe_weighted_mean(value, weight),
      .groups = "drop"
    )

  categorical_summary <- balance_data %>%
    mutate(imd_missing = if_else(is.na(imd_decile), "Missing", "Observed")) %>%
    pivot_longer(c(gender, ethnicity_5cat, imd_missing),
                 names_to = "variable", values_to = "category") %>%
    group_by(variable, category, group) %>%
    summarise(n = n(), weighted_n = sum(weight), .groups = "drop") %>%
    complete(nesting(variable, category), group = c("Cases", "Controls"),
             fill = list(n = 0L, weighted_n = 0)) %>%
    group_by(variable, group) %>%
    mutate(percent = 100 * n / sum(n),
           weighted_percent = 100 * weighted_n / sum(weighted_n)) %>%
    ungroup()

  # Fixed case-group SD denominator for gross and weighted SMDs.
  # These are the retained matched cases, not the original pre-match population.
  # IMD mean/SD/SMD are conditional on IMD being observed; missingness is separate.
  numeric_smd <- numeric_summary %>%
    select(variable, group, mean, weighted_mean, sd) %>%
    pivot_wider(names_from = group, values_from = c(mean, weighted_mean, sd)) %>%
    transmute(
      variable, category = NA_character_, denominator = sd_Cases,
      unweighted_difference = mean_Cases - mean_Controls,
      weighted_difference = weighted_mean_Cases - weighted_mean_Controls
    )

  categorical_smd <- categorical_summary %>%
    select(variable, category, group, percent, weighted_percent) %>%
    pivot_wider(names_from = group, values_from = c(percent, weighted_percent)) %>%
    transmute(
      variable, category,
      denominator = sqrt((percent_Cases / 100) * (1 - percent_Cases / 100)),
      unweighted_difference = (percent_Cases - percent_Controls) / 100,
      weighted_difference = (weighted_percent_Cases - weighted_percent_Controls) / 100
    )

  smd_summary <- bind_rows(numeric_smd, categorical_smd) %>%
    mutate(
      # Zero variance makes this SMD undefined, even if both groups agree.
      unweighted_smd = if_else(denominator > 0,
                              unweighted_difference / denominator, NA_real_),
      weighted_smd = if_else(denominator > 0,
                            weighted_difference / denominator, NA_real_),
      absolute_weighted_smd = abs(weighted_smd)
    ) %>%
    select(variable, category, unweighted_difference, weighted_difference,
           unweighted_smd, weighted_smd, absolute_weighted_smd)

  # Numeric differences: years or deciles. Categorical differences: proportions
  # (multiply by 100 for percentage points). Positive means higher in cases.
  print(numeric_summary, n = Inf, width = Inf)
  print(categorical_summary, n = Inf, width = Inf)
  print(smd_summary, n = Inf, width = Inf)
}

# 3. Individual pairs: both ages refer to the case's index date.
pair_comparison <- balance_controls %>%
  left_join(
    balance_cases %>% transmute(
      matched_case_patid = patid, case_index = as.Date(index_date),
      case_dob = as.Date(dob), case_gender = as.character(gender),
      case_ethnicity = coalesce(as.character(ethnicity_5cat), "Missing"),
      case_imd = as.numeric(na_if(as.character(imd_decile), "Missing"))
    ),
    by = "matched_case_patid"
  ) %>%
  mutate(
    age_gap_years = abs(as.numeric(as.Date(dob) - case_dob)) / 365.25,
    same_gender = as.character(gender) == case_gender,
    # Two Missing ethnicity values count as recording-category agreement.
    same_ethnicity = coalesce(as.character(ethnicity_5cat), "Missing") ==
      case_ethnicity,
    imd_decile = as.numeric(na_if(as.character(imd_decile), "Missing")),
    same_imd_missingness = is.na(imd_decile) == is.na(case_imd),
    both_imd_observed = !is.na(imd_decile) & !is.na(case_imd),
    imd_gap = abs(imd_decile - case_imd)
  )

stopifnot(all(as.Date(pair_comparison$matched_case_index_date) ==
                pair_comparison$case_index))

# Pair summaries give each pair equal weight: larger sets contribute more.
# IMD gaps use only pairs with both IMD values observed.
pair_quality <- pair_comparison %>% summarise(
  n_pairs = n(),
  mean_age_gap = safe_mean(age_gap_years),
  median_age_gap = median(age_gap_years, na.rm = TRUE),
  p95_age_gap = quantile(age_gap_years, 0.95, na.rm = TRUE),
  largest_age_gap = safe_max(age_gap_years),
  percent_same_gender = 100 * safe_mean(same_gender),
  percent_same_ethnicity = 100 * safe_mean(same_ethnicity),
  percent_same_imd_missingness = 100 * safe_mean(same_imd_missingness),
  n_pairs_with_both_imd_observed = sum(both_imd_observed),
  mean_imd_gap_observed = safe_mean(imd_gap),
  largest_imd_gap_observed = safe_max(imd_gap)
)
print(pair_quality, width = Inf)

# Matched cases vs unmatched cases (is there a systematic difference in age or gender for those who could not be matched?)

balance_cases %>%
  mutate(
    match_status = if_else(n_controls > 0, "Matched", "Unmatched"),
    age = as.numeric(as.Date(index_date) - as.Date(dob)) / 365.25
  ) %>%
  group_by(match_status) %>%
  summarise(
    n = n(),
    missing_age = sum(is.na(age)),
    mean_age = mean(age, na.rm = TRUE),
    median_age = median(age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(width = Inf)


# Interpretation:
# - numeric_summary: unweighted mean/SD/median, weighted mean, missingness.
# - categorical_summary: unweighted and weighted percentages, including missing IMD.
# - smd_summary: signed differences and SMDs; absolute_weighted_smd is the balance metric.
# - pair_comparison: individual linked pairs for inspection.
# - pair_quality: aggregate pair-level age gaps, category agreement and IMD gaps.
# SMD denominator is the SD in the retained matched cases (binary SD for categories).
# Zero/undefined case SD gives NA SMD; inspect the raw differences in those rows.
# IMD numerical balance is conditional on observation; assess missingness separately.
# These checks describe only the supplied matched sample, not unmatched cases.
