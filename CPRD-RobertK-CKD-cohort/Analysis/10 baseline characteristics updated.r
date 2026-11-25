#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(ggplot2)
library(quarto)
library(flextable)
rm(list=ls())

# Load in combined cohort (NB may need to update cohort names if the files are updated)

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/18112026_combined_ckd_cohort_dm.Rda")
combined_2024 <- get(loaded_file[1])

##########

# Check for pre-index dementia (this should have been excluded prior to matching)

pre_dementia_check <- combined_2024 %>%
  group_by(is_case) %>%
  summarise(
    n_records = n(),
    n_pre_index_dementia =
      sum(!is.na(pre_index_date_earliest_alldementia)),
    .groups = "drop"
  )

print(pre_dementia_check)

combined_summary <- combined_2024 %>%
  filter(is.na(pre_index_date_earliest_alldementia) & is.na(pre_index_date_latest_dementia_medications))

# Show how many records were excluded
tibble(
  original_records = nrow(combined_2024),
  included_records = nrow(combined_summary),
  excluded_records = nrow(combined_2024) - nrow(combined_summary)
)

# Create dementia, duration and follow-up variables

combined_summary <- combined_summary %>%
  mutate(
    post_index_alldementia =
      as.integer(!is.na(post_index_date_first_alldementia) | !(is.na(post_index_date_first_dementia_medications))),

    follow_up_end_date = pmin(gp_end_date, death_date, hes_end_date, na.rm = TRUE),

    follow_up_length = as.numeric(
      difftime(follow_up_end_date, index_date, units = "days")
    ),

    follow_up_years = follow_up_length / 365.25
  )

# Prepare variables for the summary table

summary_data <- combined_summary %>%
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

    # Missing dm_at_index denotes no diabetes in the merged cohort.
    diabetes_present = coalesce(dm_at_index == 1, FALSE),
    diabetes_duration_group = case_when(
      !diabetes_present ~ "No diabetes",
      TRUE ~ as.character(cut(
        dm_dur_all,
        breaks = c(0, 10, 20, 30, Inf),
        labels = c("0–<10", "10–<20", "20–<30", "30+"),
        right = FALSE
      ))
    ),

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

    `Never smoker` =
      n_percent(never_smoker),

    `Hypertension` =
      n_percent(hypertension_present),

    `Diabetes n (%)` =
      n_percent(diabetes_present),

    `No diabetes` =
      n_percent(!diabetes_present),

    `Diabetes duration 0–<10 years` =
      n_percent(diabetes_duration_group %in% "0–<10"),

    `Diabetes duration 10–<20 years` =
      n_percent(diabetes_duration_group %in% "10–<20"),

    `Diabetes duration 20–<30 years` =
      n_percent(diabetes_duration_group %in% "20–<30"),

    `Diabetes duration 30+ years` =
      n_percent(diabetes_duration_group %in% "30+"),

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
    `Never smoker`,
    `Hypertension`,
    `Diabetes n (%)`,
    `No diabetes`,
    starts_with("Diabetes duration"),
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
saveRDS(summary_table, "C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/18112026_baseline_characteristics_summary.rds")

quarto::quarto_render("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Outputs/Analysis outputs.qmd")


# Report missing data separately

missing_data_table <- summary_data %>%
  group_by(cohort, dementia_group) %>%
  summarise(
    n_records = n(),
    missing_age = sum(is.na(index_date_age)),
    missing_gender = sum(is.na(male)),
    missing_ethnicity = sum(is.na(white)),
    missing_smoking = sum(is.na(never_smoker)),
    missing_or_invalid_duration_among_diabetes =
      sum(diabetes_present & is.na(diabetes_duration_group)),
    missing_imd = sum(is.na(imd_group)),
    missing_follow_up = sum(is.na(follow_up_years)),
    .groups = "drop"
  )

print(missing_data_table, n = Inf, width = Inf)