#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(ggplot2)
rm(list=ls())

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/05102026_matched_ckd_cohort_dm.Rda")

diabetes_2024 <- get(loaded_file[1])

# Keep individuals aged 18 or over
diabetes_2024 <- diabetes_2024 %>%
  filter(index_date_age >= 18)

dim(diabetes_2024)
names(diabetes_2024)
head(diabetes_2024)
View(diabetes_2024)
summary(diabetes_2024)

# Column names containing "ckd", regardless of capitalisation
names(diabetes_2024)[
  grepl("ckd", names(diabetes_2024), ignore.case = TRUE)
]

diabetes_2024 %>%
  select(patid, index_date_age, gender, dob, pracid, dm_dur_all, index_date_ckd_dur_all, preckdstage, pre_index_date_earliest_ckd5_code, post_index_date_first_ckd5_code) %>% 
  head(50)

# Number of missing values in each column, largest first
sort(colSums(is.na(diabetes_2024)), decreasing = TRUE)

# Code for summarising variables
summary(diabetes_2024$index_date_age)
summary(diabetes_2024$index_date)
summary(diabetes_2024$dm_dur_all)

diabetes_2024 %>%
  summarise(
    n_observed = sum(!is.na(index_date)),
    n_missing = sum(is.na(index_date)),
    mean_age = mean(index_date, na.rm = TRUE),
    sd_age = sd(index_date, na.rm = TRUE),
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
    summarise(index_date_age = mean(index_date_age, na.rm = TRUE)) %>%
    #mutate(percent = round(100 * n / sum(n), 1))

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