# Balance and pair-quality summaries only
# Run AFTER your matching script, in the same R session.
# Inputs: local matched_cases and matched_controls tables.
# No matching, database writes or changes to these input tables are performed.
# Replace your existing "Checking quality of matches" section with this script,
# or source this file after the matching/output sections finish.
# No cobalt dependency. Missing IMD remains eligible.
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

# Each matched case contributes weight 1; its controls together contribute 1.
# Unmatched cases are counted above but excluded from balance comparisons.
balance_data <- bind_rows(
  balance_cases %>%
    filter(n_controls > 0) %>%
    transmute(group = "Cases", dob, gender, ethnicity_5cat, imd_decile,
              index_date, weight = 1),
  balance_controls %>%
    left_join(set_sizes, by = "matched_case_patid") %>%
    transmute(group = "Controls", dob, gender, ethnicity_5cat, imd_decile,
              index_date = matched_case_index_date, weight = 1 / n_controls)
) %>%
  mutate(
    age = as.numeric(as.Date(index_date) - as.Date(dob)) / 365.25,
    gender = as.character(gender),
    ethnicity_5cat = coalesce(as.character(ethnicity_5cat), "Missing"),
    imd_decile = as.numeric(na_if(as.character(imd_decile), "Missing"))
  ) %>%
  select(group, age, gender, ethnicity_5cat, imd_decile, weight)

# IMD alone may be NA; the 5.5 matching placeholder is never summarised.
stopifnot(
  all(complete.cases(select(balance_data, -imd_decile))),
  all(is.finite(balance_data$weight)), all(balance_data$weight > 0)
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
