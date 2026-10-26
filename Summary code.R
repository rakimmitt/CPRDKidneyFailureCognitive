#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(ggplot2)
library(cmprsk)
rm(list=ls())

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/20102026_matched_ckd_cohort_dm.Rda")

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

diabetes_2024 %>%
filter(is_case == 0) %>%
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