#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(ggplot2)
rm(list=ls())

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/05102026_matched_ckd_cohort_dm.Rda")

diabetes_2024 <- get(loaded_file[1])

# Identify the patient IDs of cases aged under 18
patid_under_18 <- diabetes_2024 %>%
  filter(is_case == 1, index_date_age < 18) %>%
  pull(patid)

# Remove under-18 individuals and controls matched to under-18 cases
diabetes_2024 <- diabetes_2024 %>%
  filter(
    index_date_age >= 18,
    !(is_case == 0 & matched_case_patid %in% patid_under_18)
  )

# Keep individuals aged 18 or over - this will remove any controls <18 matched to cases >18
diabetes_2024 <- diabetes_2024 %>%
  filter(index_date_age >= 18)


# Flag individuals with either pre-index dementia date recorded
dementia_check <- diabetes_2024 %>%
  mutate(
    pre_index_alldementia =
      !is.na(pre_index_date_latest_alldementia) |
      !is.na(pre_index_date_earliest_alldementia)
  )

# Count exclusions separately for cases and controls
exclusion_summary <- dementia_check %>%
  group_by(is_case) %>%
  summarise(
    n_original = n(),
    n_excluded = sum(pre_index_alldementia),
    n_retained = sum(!pre_index_alldementia),
    percent_excluded = 100 * n_excluded / n_original,
    .groups = "drop"
  )

exclusion_summary %>%
  print(width = Inf)

# Overall number excluded
sum(dementia_check$pre_index_alldementia)

# Exclude pre-index dementia and create variable to flag post-index dementia diagnoses
dementia_2024 <- dementia_check %>%
  filter(!pre_index_alldementia) %>%
  mutate(
    post_index_alldementia = as.integer(
        !is.na(post_index_date_first_alldementia)
      ),
    time_to_alldementia = as.numeric(
      difftime(
        post_index_date_first_alldementia,
        index_date,
        units = "days"
      )
    )
  )

# Count post-index dementia diagnoses by case/control group
dementia_2024 %>%
  group_by(is_case) %>%
  summarise(
    n = n(),
    n_post_index_alldementia = sum(post_index_alldementia),
    percent_post_index_alldementia =
      100 * n_post_index_alldementia / n,
    .groups = "drop"
  ) %>%
  print(width = Inf)

# Summarise time to dementia among those diagnosed
dementia_2024 %>%
  filter(post_index_alldementia == 1) %>%
  group_by(is_case) %>%
  summarise(
    n_diagnosed = n(),
    n_with_time = sum(!is.na(time_to_alldementia)),
    mean_days = mean(time_to_alldementia, na.rm = TRUE),
    sd_days = sd(time_to_alldementia, na.rm = TRUE),
    median_days = median(time_to_alldementia, na.rm = TRUE),
    q1_days = quantile(time_to_alldementia, 0.25, na.rm = TRUE),
    q3_days = quantile(time_to_alldementia, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# Flag for pre-index hypertension
  dementia_2024 <- dementia_2024 %>%
  mutate(
    pre_index_hypertension = !is.na(pre_index_date_latest_alldementia)
  )

### Baseline characteristics table ###

  # Number and percentage among those with non-missing information
n_percent <- function(condition) {
  denominator <- sum(!is.na(condition))

  if (denominator == 0) return(NA_character_)

  sprintf(
    "%s (%.1f%%)",
    format(sum(condition, na.rm = TRUE),
           big.mark = ",", trim = TRUE),
    100 * sum(condition, na.rm = TRUE) / denominator
  )
}

# Mean (SD)
mean_sd <- function(x) {
  if (all(is.na(x))) return(NA_character_)

  sprintf(
    "%.1f (%.1f)",
    mean(x, na.rm = TRUE),
    sd(x, na.rm = TRUE)
  )
}

# Median [lower quartile, upper quartile]
median_iqr <- function(x) {
  if (all(is.na(x))) return(NA_character_)

  sprintf(
    "%.1f [%.1f, %.1f]",
    median(x, na.rm = TRUE),
    quantile(x, 0.25, na.rm = TRUE),
    quantile(x, 0.75, na.rm = TRUE)
  )
}

# Prepare variables 
# NB smoker has been silenced as currently missing data

table_data <- dementia_2024 %>%
  mutate(
    male = gender == 1,
    white = ethnicity_5cat == 0,
    diabetes_duration_years = dm_dur_all,
    #never_smoker = smoking_cat == "Never smoker",
    baseline_hypertension = pre_index_hypertension == 1,

    group = factor(
      case_when(
        is_case == 1 & post_index_alldementia == 1 ~
          "aCKD: post-index dementia",
        is_case == 1 & post_index_alldementia == 0 ~
          "aCKD: no dementia",
        is_case == 0 & post_index_alldementia == 1 ~
          "Control: post-index dementia",
        is_case == 0 & post_index_alldementia == 0 ~
          "Control: no dementia"
      ),
      levels = c(
        "aCKD: post-index dementia",
        "aCKD: no dementia",
        "Control: post-index dementia",
        "Control: no dementia"
      )
    )
  )

# Calculate the table ----------------------------------------------

summary_table <- table_data %>%
  filter(!is.na(group)) %>%
  group_by(group, .drop = FALSE) %>%
  summarise(
    `n` = format(n(), big.mark = ",", trim = TRUE),

    `Mean age (SD)` = mean_sd(index_date_age),

    `Male, n (%)` = n_percent(male),

    `White ethnicity, n (%)` = n_percent(white),

    `Diabetes duration, years: median [Q1, Q3]` =
      median_iqr(diabetes_duration_years),

    #`Never smoker, n (%)` = n_percent(never_smoker),

    `Hypertension, n (%)` = n_percent(baseline_hypertension),

    `IMD 1–2 (least deprived), n (%)` =
      n_percent(imd_decile >= 1 & imd_decile <= 2),

    `IMD 3–4, n (%)` =
      n_percent(imd_decile >= 3 & imd_decile <= 4),

    `IMD 5–6, n (%)` =
      n_percent(imd_decile >= 5 & imd_decile <= 6),

    `IMD 7–8, n (%)` =
      n_percent(imd_decile >= 7 & imd_decile <= 8),

    `IMD 9–10 (most deprived), n (%)` =
      n_percent(imd_decile >= 9 & imd_decile <= 10),

    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -group,
    names_to = "Characteristic",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = value
  )

View(summary_table)

summary_table %>%
  print(n = Inf, width = Inf)