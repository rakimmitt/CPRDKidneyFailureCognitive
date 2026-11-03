#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(ggplot2)
library(cmprsk)
rm(list=ls())

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/03112026_matched_ckd_cohort_dm.Rda")

diabetes_2024 <- get(loaded_file[1])

dim(diabetes_2024)
names(diabetes_2024)
head(diabetes_2024)
View(diabetes_2024)
summary(diabetes_2024)

diabetes_2024 %>%
  select(patid, index_date_age, gender, dob, pracid, dm_dur_all, index_date_ckd_dur_all, preckdstage, pre_index_date_earliest_ckd5_code, post_index_date_first_ckd5_code) %>% 
  head(50)

# Number of missing values in each column, largest first
sort(colSums(is.na(diabetes_2024)), decreasing = TRUE)

# Code for summarising variables
summary(diabetes_2024$index_date_age)
summary(diabetes_2024$index_date)
summary(diabetes_2024$dm_dur_all)
summary(diabetes_2024$preegfr)

diabetes_2024 %>%
filter(is_case ==0) %>%
summary(preegfr)

diabetes_2024 <- diabetes_2024 %>%
  mutate(
    registration_duration_days = as.numeric(
      difftime(gp_end_date, regstartdate, units = "days")
    )
  )

# Create a follow-up end date and look at follow-up duration (in years)
diabetes_2024 <- diabetes_2024 %>%
  mutate(
    follow_up_end_date = pmin(
      gp_end_date, death_date, hes_end_date,
      na.rm = TRUE
    ),
    follow_up_length = as.numeric(
      difftime(follow_up_end_date, index_date, units = "days")
    ),
    follow_up_years = follow_up_length / 365.25
  )
summary(diabetes_2024$follow_up_years)

diabetes_2024 %>%
  summarise(
    n_observed = sum(!is.na(index_date)),
    n_missing = sum(is.na(index_date)),
    mean_age = mean(index_date_age, na.rm = TRUE),
    sd_age = sd(index_date_age, na.rm = TRUE),
    median_age = median(index_date, na.rm = TRUE),
    iqr_age = IQR(index_date, na.rm = TRUE)
  )

  diabetes_2024 %>%
    summarise(
      across(
        c(index_date_age, preegfr, prebmi),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE),
          missing = ~ sum(is.na(.x))
        )
      )
    )
  
  # Counts, including missing values
  table(diabetes_2024$gender, useNA = "ifany")
  table(diabetes_2024$index_date_ckd_stage, useNA = "ifany")
  table(diabetes_2024$is_case, useNA = "ifany")

# Counts with a condition

  diabetes_2024 %>%
    filter(is_case == 0) %>%
    summarise(preegfr = mean(preegfr, na.rm = TRUE)) %>%
    #mutate(percent = round(100 * n / sum(n), 1))

column_name <- "preegfr"
diabetes_summary %>%
filter(is_case == 1) %>%
  summarise(
    n_observed = sum(!is.na(.data[[column_name]])),
    n_missing = sum(is.na(.data[[column_name]])),
    mean = mean(.data[[column_name]], na.rm = TRUE),
    sd = sd(.data[[column_name]], na.rm = TRUE),
    median = median(.data[[column_name]], na.rm = TRUE),
    iqr = IQR(.data[[column_name]], na.rm = TRUE),
    max = max(.data[[column_name]], na.rm = TRUE),
    min = min(.data[[column_name]], na.rm = TRUE)
  )

  diabetes_2024 %>%
    filter(is_case == 0) %>%
    summarise(dm_dur_all = mean(dm_dur_all, na.rm = TRUE)) %>%

  
  # Counts and percentages in a single table
  diabetes_2024 %>%
    count(gender, name = "n") %>%
    mutate(percent = round(100 * n / sum(n), 1))
  
  # Counts for each combination
  table(
    diabetes_2024$index_date_ckd_stage,
    diabetes_2024$post_index_date_first_dementia,
    useNA = "ifany"
  )

    table(
    diabetes_2024$is_case,
    diabetes_2024$dm_dur_all,
    useNA = "ifany"
  )
  
  
  diabetes_2024 %>%
    group_by(preckdstage) %>%
    summarise(
      n = n(),
      age_missing = sum(is.na(index_date_age)),
      mean_age = mean(index_date_age, na.rm = TRUE),
      sd_age = sd(index_date_age, na.rm = TRUE),
      median_age = median(index_date_age, na.rm = TRUE),
      .groups = "drop"
    )

      diabetes_2024 %>%
    group_by(is_case) %>%
    summarise(
      n = n(),
      age_missing = sum(is.na(index_date_age)),
      mean_age = mean(index_date_age, na.rm = TRUE),
      sd_age = sd(index_date_age, na.rm = TRUE),
      median_age = median(index_date_age, na.rm = TRUE),
      mean_dm_dur_all = mean(dm_dur_all, na.rm = TRUE),
      sd_dm_dur_all = sd(dm_dur_all, na.rm = TRUE),
      median_post_index_date_first_dementia = median(post_index_date_first_dementia, na.rm = TRUE),
      .groups = "drop",
    ) %>% print(width = Inf, length = Inf)

    diabetes_2024 %>%
    filter(is.na(pre_index_date_earliest_alldementia)) %>%
    group_by(is_case) %>%
    summarise(
      n = n(),
      dementia_diagnosis = sum(!is.na(post_index_date_first_alldementia)),
      dementia_percent = 100 * dementia_diagnosis / n,
      mean_dm_dur_all = mean(dm_dur_all, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(diabetes_2024, aes(x = index_date_age)) +
    geom_histogram(binwidth = 5, fill = "steelblue", colour = "white") +
    labs(x = "Age (years)", y = "Number of records") +
    theme_minimal()
  
  ggplot(diabetes_2024, aes(x = factor(preckdstage), y = index_date_age)) +
    geom_boxplot() +
    labs(x = "CKD stage", y = "Age (years)") +
    theme_minimal()









diabetes_summary %>%
filter(is_case ==0) %>%
summary(preegfr)


low_egfr_controls <- diabetes_2024 %>%
  filter(is_case == 0, !is.na(preegfr), preegfr < 60) %>%
  select(
    any_of(c(
      "match_record_id", "patid", "matched_case_patid",
      "index_date", "preegfr", "preegfrdate", "preegfrdatediff",
      "index_date_ckd_stage",
      "ckd_stage_3_start_date", "advanced_ckd_index_date"
    ))
  ) %>%
  arrange(preegfr)

View(low_egfr_controls)



continuous_vars <- c("age_at_index", "dm_dur_all", "prebmi")
categorical_vars <- c("gender", "smoking_cat", "imd_decile", "ethnicity_5cat")

# Separate histograms for continuous variables
for (v in continuous_vars) {
  p <- ggplot(diabetes_2024, aes(x = .data[[v]])) +
    geom_histogram(
      bins = 30,
      fill = "steelblue",
      colour = "white",
      na.rm = TRUE
    ) +
    labs(
      title = v,
      subtitle = paste("Missing:", sum(is.na(diabetes_2024[[v]]))),
      x = v,
      y = "Number of patients"
    ) +
    theme_minimal()

  print(p)
}

# Separate bar charts for categorical variables
# Missing values appear as an NA category
for (v in categorical_vars) {
  p <- ggplot(diabetes_2024, aes(x = factor(.data[[v]]))) +
    geom_bar(fill = "steelblue") +
    scale_x_discrete(na.translate = TRUE) +
    labs(
      title = v,
      x = v,
      y = "Number of patients"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p)
}




# Dementia medication concordance

dementia_check <- diabetes_2024 %>%
  transmute(
    diagnosis_date = as.Date(post_index_date_first_alldementia),
    medication_date = as.Date(post_index_date_first_dementia_medications),

    diagnosis_present = !is.na(diagnosis_date),
    medication_present = !is.na(medication_date),

    # Positive = medication recorded AFTER diagnosis
    # Negative = medication recorded BEFORE diagnosis
    gap_days = as.numeric(medication_date - diagnosis_date),

    concordance_group = case_when(
      diagnosis_present & medication_present  ~ "Both",
      diagnosis_present & !medication_present ~ "Diagnosis only",
      !diagnosis_present & medication_present ~ "Medication only",
      TRUE                                   ~ "Neither"
    )
  )


# 2 × 2 table
with(
  dementia_check,
  table(
    Diagnosis = factor(
      diagnosis_present,
      levels = c(FALSE, TRUE),
      labels = c("Absent", "Present")
    ),
    Medication = factor(
      medication_present,
      levels = c(FALSE, TRUE),
      labels = c("Absent", "Present")
    )
  )
)

# Counts and percentages in each group
dementia_check %>%
  count(concordance_group, name = "n") %>%
  mutate(percent = round(100 * n / sum(n), 1))

# Returns NA if the denominator is zero
pct <- function(numerator, denominator) {
  if (denominator == 0) return(NA_real_)
  round(100 * numerator / denominator, 1)
}

dementia_check %>%
  summarise(
    n_patients = n(),
    n_diagnosis = sum(diagnosis_present),
    n_medication = sum(medication_present),
    n_both = sum(diagnosis_present & medication_present),
    n_either = sum(diagnosis_present | medication_present),

    overall_agreement_pct = pct(
      sum(diagnosis_present == medication_present), n()
    ),

    # Among patients with a diagnosis, how many have medication?
    medication_among_diagnosed_pct = pct(n_both, n_diagnosis),

    # Among patients with medication, how many have a diagnosis?
    diagnosis_among_medicated_pct = pct(n_both, n_medication),

    # Overlap among patients with either marker
    both_among_either_pct = pct(n_both, n_either)
  )

  both_dates <- dementia_check %>%
  filter(diagnosis_present & medication_present)

# Numerical summary; only run if there are paired dates
if (nrow(both_dates) > 0) {

  gap_summary <- both_dates %>%
    summarise(
      n = n(),
      median_gap_days = median(gap_days),
      q1_gap_days = unname(quantile(gap_days, 0.25)),
      q3_gap_days = unname(quantile(gap_days, 0.75)),
      min_gap_days = min(gap_days),
      max_gap_days = max(gap_days),
      median_absolute_gap_days = median(abs(gap_days)),
      within_30_days_pct = pct(sum(abs(gap_days) <= 30), n()),
      within_90_days_pct = pct(sum(abs(gap_days) <= 90), n()),
      within_365_days_pct = pct(sum(abs(gap_days) <= 365), n())
    )

  print(gap_summary)

  # Direction of the gap
  timing_summary <- both_dates %>%
    mutate(
      timing = case_when(
        gap_days < 0 ~ "Medication before diagnosis",
        gap_days == 0 ~ "Same day",
        gap_days > 0 ~ "Medication after diagnosis"
      )
    ) %>%
    count(timing, name = "n") %>%
    mutate(percent = round(100 * n / sum(n), 1))

  print(timing_summary)

  # Distribution of signed time gaps
  p <- ggplot(both_dates, aes(x = gap_days)) +
    geom_histogram(
      bins = 50,
      fill = "steelblue",
      colour = "white"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Time between first recorded dementia diagnosis and medication",
      subtitle = "Patients with both post-index dates",
      x = "Days: medication date minus diagnosis date",
      y = "Number of patients"
    ) +
    theme_minimal()

  print(p)

} else {
  message("No patients have both dates.")
}
