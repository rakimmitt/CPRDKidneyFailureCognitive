#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(ggplot2)
rm(list=ls())

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/28102026_matched_ckd_cohort_dm.Rda")

diabetes_2024 <- get(loaded_file[1])

##########

# Check for pre-index dementia (this should have been excluded prior to matching)

pre_dementia_check <- diabetes_2024 %>%
  group_by(is_case) %>%
  summarise(
    n_records = n(),
    n_pre_index_dementia =
      sum(!is.na(pre_index_date_earliest_alldementia)),
    .groups = "drop"
  )

print(pre_dementia_check)

# Restrict to participants with diabetes at the index date (and keep diabetes_2024 unchanged)

diabetes_summary <- diabetes_2024 %>%
  filter(dm_at_index == 1, dm_dur_all > 0)

# Show how many records were excluded
tibble(
  original_records = nrow(diabetes_2024),
  included_records = nrow(diabetes_summary),
  excluded_records = nrow(diabetes_2024) - nrow(diabetes_summary)
)

# Create dementia, duration and follow-up variables

diabetes_summary <- diabetes_summary %>%
  mutate(
    post_index_alldementia =
      as.integer(!is.na(post_index_date_first_alldementia)),

    follow_up_end_date = pmin(gp_end_date, death_date, hes_end_date, na.rm = TRUE),

    follow_up_length = as.numeric(
      difftime(follow_up_end_date, index_date, units = "days")
    ),

    follow_up_years = follow_up_length / 365.25
  )

# Prepare variables for the summary table

summary_data <- diabetes_summary %>%
  mutate(
    cohort = case_when(
      is_case == 1 ~ "Advanced CKD",
      is_case == 0 ~ "Control",
      TRUE ~ NA_character_
    ),

    dementia_group = if_else(
      post_index_alldementia == 1,
      "Post-index dementia",
      "No dementia"
    ),

    male = as.character(gender) == 1,

# Treat the ethnicity "Missing" category as missing data
    white = na_if(
      as.character(ethnicity_5cat), "Missing"
    ) == 0,

# Use QRISK2 because smoking_cat has no explicit never-smoker category
    never_smoker = qrisk2_smoking_cat == 0,

# Baseline hypertension: at least one pre-index diagnosis
    hypertension_present =
      !is.na(pre_index_date_earliest_hypertension),

    imd_numeric = as.numeric(as.character(imd_decile)),

    imd_group = case_when(
      imd_numeric %in% 1:2  ~ "1–2",
      imd_numeric %in% 3:4  ~ "3–4",
      imd_numeric %in% 5:6  ~ "5–6",
      imd_numeric %in% 7:8  ~ "7–8",
      imd_numeric %in% 9:10 ~ "9–10",
      TRUE ~ NA_character_
    ),

    column = paste(cohort, dementia_group, sep = ": ")
  )

if (any(is.na(summary_data$cohort))) {
  stop("Some records have missing or unexpected is_case values.")
}

# Formatting functions

# Count (%) using all records in the subgroup as denominator

n_percent <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_character_)
  }

  count <- sum(x, na.rm = TRUE)

  sprintf(
    "%s (%.1f%%)",
    format(count, big.mark = ",", trim = TRUE),
    100 * count / length(x)
  )
}

# Mean (SD), excluding missing values
mean_sd <- function(x) {
  if (all(is.na(x))) {
    return(NA_character_)
  }

  sprintf(
    "%.1f (%.1f)",
    mean(x, na.rm = TRUE),
    sd(x, na.rm = TRUE)
  )
}

# Median [lower quartile, upper quartile], excluding missing values
median_iqr <- function(x) {
  if (all(is.na(x))) {
    return(NA_character_)
  }

  q <- quantile(
    x,
    probs = c(0.25, 0.5, 0.75),
    na.rm = TRUE
  )

  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

# Calculate cohort totals 

cohort_totals <- summary_data %>%
  group_by(cohort) %>%
  summarise(
    cohort_n = n(),
    n_unique_participants = n_distinct(patid),
    .groups = "drop"
  )

print(cohort_totals)


# Calculate statistics for the four groups

group_statistics <- summary_data %>%
  group_by(cohort, dementia_group, column) %>%
  summarise(
    subgroup_n = n(),

    `Mean age (SD)` =
      mean_sd(index_date_age),

    `Male` =
      n_percent(male),

    `White ethnicity` =
      n_percent(white),

    `Median diabetes duration, years [IQR]` =
      median_iqr(dm_dur_all),

    `Never smoker` =
      n_percent(never_smoker),

    `Hypertension` =
      n_percent(hypertension_present),

    `IMD decile 1–2` =
      n_percent(imd_group == "1–2"),

    `IMD decile 3–4` =
      n_percent(imd_group == "3–4"),

    `IMD decile 5–6` =
      n_percent(imd_group == "5–6"),

    `IMD decile 7–8` =
      n_percent(imd_group == "7–8"),

    `IMD decile 9–10` =
      n_percent(imd_group == "9–10"),

    `Median follow-up duration, years [IQR]` =
      median_iqr(follow_up_years),

    .groups = "drop"
  ) %>%
  left_join(
    cohort_totals %>% select(cohort, cohort_n),
    by = "cohort"
  ) %>%
  mutate(
    # For this row, the percentage denominator is the whole cohort.
    n = sprintf(
      "%s (%.1f%%)",
      format(subgroup_n, big.mark = ",", trim = TRUE),
      100 * subgroup_n / cohort_n
    )
  )

# Arrange the final table

column_order <- c(
  "Advanced CKD: Post-index dementia",
  "Advanced CKD: No dementia",
  "Control: Post-index dementia",
  "Control: No dementia"
)

summary_table <- group_statistics %>%
  select(
    column,
    n,
    `Mean age (SD)`,
    `Male`,
    `White ethnicity`,
    `Median diabetes duration, years [IQR]`,
    `Never smoker`,
    `Hypertension`,
    starts_with("IMD decile"),
    `Median follow-up duration, years [IQR]`
  ) %>%
  mutate(
    column = factor(column, levels = column_order)
  ) %>%
  pivot_longer(
    cols = -column,
    names_to = "Characteristic",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = column,
    values_from = Value,
    names_expand = TRUE
  ) %>%
  select(Characteristic, all_of(column_order))

print(summary_table, n = Inf, width = Inf)

View(summary_table)

# Report missing data separately

missing_data_table <- summary_data %>%
  group_by(cohort, dementia_group) %>%
  summarise(
    n_records = n(),
    missing_age = sum(is.na(index_date_age)),
    missing_gender = sum(is.na(male)),
    missing_ethnicity = sum(is.na(white)),
    missing_diabetes_duration = sum(is.na(dm_dur_all)),
    missing_smoking = sum(is.na(never_smoker)),
    missing_imd = sum(is.na(imd_group)),
    missing_follow_up = sum(is.na(follow_up_years)),
    .groups = "drop"
  )

print(missing_data_table, n = Inf, width = Inf)