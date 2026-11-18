#Setup
library(tidyverse)
library(dplyr)
library(aurum)
library(EHRBiomarkr)
library(ggplot2)

# Load in both cohorts - NB may need to update cohort names if the files are updated

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/18112026_matched_ckd_cohort_dm.Rda")
diabetes_2024 <- get(loaded_file[1])

loaded_file <- load("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/12112026_matched_ckd_cohort_nondm.Rda")
nondiabetes_2024 <- get(loaded_file[1])

#############

# Check non-overlapping columns between the two datasets

column_inventory <- function(data, type_column) {
  result <- tibble(
    column = names(data),
    type = vapply(
      data,
      function(x) paste(class(x), collapse = "/"),
      character(1)
    )
  )

  names(result)[2] <- type_column
  result
}

column_comparison <- full_join(
  column_inventory(diabetes_2024, "type_diabetes"),
  column_inventory(nondiabetes_2024, "type_nondiabetes"),
  by = "column"
) %>%
  mutate(
    location = case_when(
      is.na(type_diabetes) ~ "Non-diabetes dataset only",
      is.na(type_nondiabetes) ~ "Diabetes dataset only",
      TRUE ~ "Both"
    ),
    same_type = type_diabetes == type_nondiabetes
  ) %>%
  arrange(location, column)

print(column_comparison, n = Inf)
View(column_comparison)

# Specifically inspect columns that are absent or have different types
column_comparison %>%
  filter(location != "Both" | !same_type) %>%
  print(n = Inf)

###### Prepare working copies for merging

dm <- diabetes_2024
nondm <- nondiabetes_2024

# Explicitly list wanted columns prior to merging

keep_columns <- c(
  "patid",
  "dob",
  "gender",
  "index_date",
  "index_date_age",
  "dm_at_index",
  "dm_dur_all",
  "dm_diag_date_all",
  "diabetes_diagnosis_date",
  "advanced_ckd_index_date",
  "alcohol_cat",
  "ckd_stage_3_start_date",
  "death_date",
  "earliest_all_dementia",
  "efi_n_deficits",
  "ethnicity_5cat",
  "gp_end_date",
  "height",
  "hes_end_date",
  "imd_decile",
  "index_date_ckd_stage",
  "is_case",
  "match_record_id",
  "matched_case_index_date",
  "matched_case_patid",
  "n_controls",
  "pracid",
  "preacr",
  "preacrdatediff",
  "prebmi",
  "prebmidate",
  "prebmidatediff",
  "preegfr",
  "preegfrdate",
  "preegfrdatediff",
  "prehaemoglobin",
  "prehaemoglobindate",
  "prehaemoglobindatediff",
  "prehba1c",
  "prehba1cdate",
  "prehba1cdatediff",
  "prehdl",
  "prehdldate",
  "prehdldatediff",
  "presbp",
  "presbpdate",
  "presbpdatediff",
  "pretotalcholesterol",
  "pretotalcholesteroldate",
  "pretotalcholesteroldatediff",
  "primary_death_cause1",
  "primary_death_cause2",
  "primary_death_cause3",
  "qrisk2_smoking_cat",
  "qrisk2_smoking_cat_uncoded",
  "regstartdate",
  "smoking_cat",
  "tds_2011",
  "weight",
  "with_hes",
  "pre_index_date_earliest_frailty_mild",
  "pre_index_date_earliest_frailty_moderate",
  "pre_index_date_earliest_frailty_severe",
  "pre_index_date_earliest_insulin",
  "pre_index_date_frailty_mild",
  "pre_index_date_frailty_moderate",
  "pre_index_date_frailty_severe",
  "pre_index_date_latest_insulin",
  "any_gn_cause",
  "any_nonstandard_ckd_cause",
  "any_nonstandard_ckd_cause_excl_igan",
  "any_pkd_cause",
  "pre_index_date_efi_hearing_impairment",
  "pre_index_date_efi_visual_impairment"
)

conditions <- c(
  "acutesinusitis",
  "af",
  "alldementia",
  "alzheimers",
  "anxiety_disorders",
  "boneinfection",
  "cellulitis",
  "ckd5_code",
  "ckd5_nokrt",
  "covid",
  "delirium",
  "dementia",
  "dka",
  "endocarditis",
  "eyeinfection",
  "falls",
  "genital_infection_nonspec",
  "giinfection",
  "haemodialysis",
  "hypertension",
  "incident_stroke",
  "infectiveotitisexterna",
  "influenza",
  "jointinfection",
  "lrti",
  "mci",
  "meningitis",
  "otherfungalinfection",
  "otherskininfection",
  "peritoneal_dialysis",
  "pneumonia",
  "primary_incident_stroke",
  "renalaccessinfection",
  "respiratorytractinfection",
  "sepsis",
  "skininfection",
  "solidorgantransplant",
  "stroke",
  "surgicalsiteinfection",
  "tia",
  "transplant",
  "tuberculosis",
  "urti",
  "uti",
  "vascular_dementia"
  )

condition_columns <- c(
  paste0("post_index_date_first_", conditions),
  paste0("pre_index_date_", conditions),
  paste0("pre_index_date_earliest_", conditions),
  paste0("pre_index_date_latest_", conditions)
)

medications <- c(
  "ace_inhibitors",
  "arb",
  "beta_blockers",
  "ksparing_diuretics",
  "loop_diuretics",
  "ca_channel_blockers",
  "thiazide_diuretics",
  "statins",
  "dementia_medications"
)

medication_columns <- c(
  paste0("post_index_date_first_", medications),
  paste0("pre_index_date_earliest_", medications),
  paste0("pre_index_date_latest_", medications)
)

kidney_cause <- c(
  "aav",
  "adpkd",
  "alport",
  "antigbm",
  "fabry",
  "fsgs",
  "gn_nos",
  "igan",
  "mcd",
  "membranous",
  "mpgn",
  "other_pkd",
  "sle"
  )

kidney_cause_columns <- c(
  paste0("ckd_cause_", kidney_cause),
  paste0("ckd_cause_earliest_", kidney_cause),
  paste0("pre_index_date_", kidney_cause),
  paste0("pre_index_date_earliest_", kidney_cause),
  paste0("pre_index_date_latest_", kidney_cause)
)

keep_columns <- unique(c(keep_columns, condition_columns, medication_columns, kidney_cause_columns))

# Check for requested columns that are absent from BOTH datasets.

missing_columns <- setdiff(keep_columns, union(names(dm), names(nondm)))

if (length(missing_columns) > 0L) {
  stop(
    "These requested columns are missing from both datasets: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Merge 

combined_2024 <- bind_rows(
  diabetes = dm %>% select(any_of(keep_columns)),
  nondiabetes = nondm %>% select(any_of(keep_columns)),
  .id = "source_dataset"
) %>%
  select(source_dataset, all_of(keep_columns))

# Confirm that all rows were retained.
stopifnot(
  nrow(combined_2024) == nrow(dm) + nrow(nondm)
)

combined_2024 %>%
  count(source_dataset, dm_at_index) %>%
  print(n = Inf)

# Create ordinal categorial variable for diabetes duration

duration_breaks <- c(seq(0, 30, by = 2), Inf)
duration_labels <- c(
  paste0(seq(0, 28, by = 2), "-<", seq(2, 30, by = 2), " years"),
  "30+ years"
)

combined_2024 <- combined_2024 %>%
  mutate(
    dm_at_index = coalesce(dm_at_index, 0L),
    diabetes_duration_cat = cut(
      dm_dur_all,
      breaks = duration_breaks,
      labels = duration_labels,
      right = FALSE,
      ordered_result = TRUE
    ),

    diabetes_duration_cat = case_when(
      dm_at_index == 0L ~ "No diabetes",
      dm_at_index == 1L ~ as.character(diabetes_duration_cat),
      TRUE ~ NA_character_
    ),

    diabetes_duration_cat = factor(
      diabetes_duration_cat,
      levels = c("No diabetes", duration_labels),
      ordered = TRUE
    )
  )

# Check baseline status by original dataset.
combined_2024 %>%
  count(source_dataset, dm_at_index) %>%
  print(n = Inf)

# Check category counts and missing duration.
combined_2024 %>%
  count(dm_at_index, diabetes_duration_cat, .drop = FALSE) %>%
  print(n = Inf)

combined_2024 %>%
  summarise(
    n_records = n(),
    n_unknown_diabetes_status = sum(is.na(dm_at_index)),
    n_diabetes_missing_duration = sum(
      dm_at_index == 1L & is.na(dm_dur_all),
      na.rm = TRUE
    )
  ) %>%
  print()

today <- format(Sys.Date(), "%d%m%Y")

save(combined_2024, file = paste0("C:/Users/rk535/OneDrive - University of Exeter/","CPRD/2024/Raw data/", today, "_combined_ckd_cohort_dm.Rda"))