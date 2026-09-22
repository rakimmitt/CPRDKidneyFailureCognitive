################################################################################
# ADVANCED CKD COHORT: KIDNEY ELIGIBILITY ONLY
#
# A patient qualifies through EITHER:
#   1. An existing algorithm-derived stage_4 or stage_5 date; OR
#   2. A valid occurrence from the CUSTOM ckd5 codelist.
#
# index_date = earliest available date from these three routes.
# The algorithm's dates are used as supplied; chronicity is not recalculated.
# A custom ckd5 code can qualify someone without an algorithm stage date.
# No separate ckd5_nokrt, dialysis, transplant or standard ckd5_code lists are used.
#
# Includes both diabetes types; retains the existing practice/gender exclusions.
# Does not apply age, registration-duration, dementia or matching restrictions.
# Therefore index_date is FIRST KIDNEY ELIGIBILITY, not necessarily the eventual
# analysis-entry date after other study eligibility criteria have been applied.
#
# All large patient-data operations remain on the server. Only the small local
# codelists are uploaded. Intermediate and final results are persistently cached.
# Existing upstream caches must already exist (as in the original file 01).
################################################################################

library(tidyverse)
library(aurum)

cprd <- CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/.aurum.yaml")

codelist_root <- paste0(
  "C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/",
  "Github clone/CPRDKidneyFailureCognitive/CPRD-Codelists"
)

# Preserve the hospital-data cutoff used in the supplied file 04.
# GP occurrences instead use each patient's gp_end_date from validDateLookup.
hes_end_date <- as.Date("2023-03-31")

cache_version <- "v1"
cache_name <- function(stem) paste0("rk_elig_", cache_version, "_", stem)

# 'analysis' specifies the namespace in which aurum stores/loads each cache.
input_analysis <- cprd$analysis("all_patid")
exclusion_analysis <- cprd$analysis("diabetes_cohort")
cohort_analysis <- cprd$analysis("all")
output_analysis <- cprd$analysis("ckd")

# 2. Load existing upstream tables ---------------------------------------------
ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  input_analysis$cached("ckd_stages_from_algorithm")

practice_exclusion_ids <- practice_exclusion_ids %>%
  exclusion_analysis$cached("practice_exclusion_ids")

gender_exclusion_ids <- gender_exclusion_ids %>%
  exclusion_analysis$cached("gender_exclusion_ids")

diabetes_cohort <- diabetes_cohort %>%
  cohort_analysis$cached("diabetes_cohort")

# Keep the earliest date for each stage if upstream contains repeated patients.
# On the database, MIN ignores NULL values; an entirely absent stage stays NULL.
# Cache BEFORE any downstream joins, so this work can be reused.
stage45_dates <- ckd_stages_from_algorithm %>%
  filter(!is.na(stage_4) | !is.na(stage_5)) %>%
  group_by(patid) %>%
  summarise(
    stage_4_date = min(as.Date(stage_4), na.rm = TRUE),
    stage_5_date = min(as.Date(stage_5), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  input_analysis$cached(
    cache_name("stage45_dates"), unique_indexes = "patid"
  )

# 3. Read ONLY the local custom ckd5 codelists -----------------------------------

read_ckd5_list <- function(folder, wanted_name, code_column) {
  directory <- file.path(codelist_root, folder)
  if (!dir.exists(directory)) stop("Codelist directory missing: ", directory)

  files <- list.files(directory, pattern = "\\.txt$", ignore.case = TRUE,
                      recursive = TRUE, full.names = TRUE)
  names_in_files <- files %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    str_to_lower() %>%
    str_remove("^exeter_medcodelist_") %>%
    str_remove("^exeter_")
  selected <- files[names_in_files == wanted_name]

  if (length(selected) == 0L) return(NULL)
  if (length(selected) > 1L) stop("Multiple codelists found for: ", wanted_name)

  codes <- readr::read_tsv(
    selected, col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, progress = FALSE
  ) %>%
    rename_with(str_to_lower) %>%
    select(all_of(code_column)) %>%
    filter(!is.na(.data[[code_column]]), .data[[code_column]] != "") %>%
    distinct()

  if (nrow(codes) == 0L) stop("Empty codelist: ", selected)
  codes
}

ckd5_medcodes <- read_ckd5_list("Medcodes", "ckd5", "medcodeid")
ckd5_icd10 <- read_ckd5_list("ICD10", "icd10_ckd5", "icd10")
ckd5_opcs4 <- read_ckd5_list("OPCS4", "opcs4_ckd5", "opcs4")

if (is.null(ckd5_medcodes) && is.null(ckd5_icd10) && is.null(ckd5_opcs4)) {
  stop("No custom ckd5 codelist was found in any coding system.")
}

# 4. Find and cache dated ckd5 code occurrences ---------------------------------
# inner_join finds matching coded records. copy=TRUE uploads the small codelist.
# Every raw result is cached so large source-table searches need not be repeated.
# code_events holds database references, not downloaded patient data.
code_events <- list()

if (!is.null(ckd5_medcodes)) {
  code_events$gp <- cprd$tables$observation %>%
    inner_join(ckd5_medcodes, by = "medcodeid", copy = TRUE) %>%
    transmute(patid, date = as.Date(obsdate), source = "gp") %>%
    input_analysis$cached(
      cache_name("raw_ckd5_gp"), indexes = c("patid", "date")
    )
}

if (!is.null(ckd5_icd10)) {
  # Match ICD prefixes, following file 04 (e.g. a parent code plus subcodes).
  code_events$icd10 <- cprd$tables$hesDiagnosisEpi %>%
    inner_join(ckd5_icd10, sql_on = "LHS.ICD LIKE CONCAT(icd10,'%')",
               copy = TRUE) %>%
    transmute(patid, date = as.Date(epistart), source = "hes") %>%
    input_analysis$cached(
      cache_name("raw_ckd5_icd10"), indexes = c("patid", "date")
    )
}

if (!is.null(ckd5_opcs4)) {
  code_events$opcs4 <- cprd$tables$hesProceduresEpi %>%
    inner_join(ckd5_opcs4, by = c("OPCS" = "opcs4"), copy = TRUE) %>%
    transmute(patid, date = as.Date(evdate), source = "hes") %>%
    input_analysis$cached(
      cache_name("raw_ckd5_opcs4"), indexes = c("patid", "date")
    )
}

# Reduce combines all available sources using union_all (append rows).
# With one source it returns that source. Duplicates cannot change MIN(date).
all_ckd5_events <- Reduce(dplyr::union_all, code_events)

# Apply the same code-date validity rules as file 04.
valid_ckd5_events <- all_ckd5_events %>%
  inner_join(cprd$tables$validDateLookup, by = "patid") %>%
  filter(
    !is.na(date), date >= min_dob,
    (source == "gp" & date <= gp_end_date) |
      (source == "hes" & date <= !!hes_end_date)
  ) %>%
  select(patid, date, source) %>%
  input_analysis$cached(
    cache_name("valid_ckd5_events"), indexes = c("patid", "date")
  )

# One row per patient: earliest qualifying CUSTOM ckd5 code in any source.
ckd5_code_dates <- valid_ckd5_events %>%
  group_by(patid) %>%
  summarise(ckd5_code_date = min(date, na.rm = TRUE), .groups = "drop") %>%
  input_analysis$cached(
    cache_name("ckd5_code_dates"), unique_indexes = "patid"
  )

# 5. Define the first kidney-eligibility date -----------------------------------
# Turn the three date columns into rows, omitting absent dates. This avoids
# placeholder dates such as 2050-01-01 and preserves patients found by codes only.
entry_dates <- stage45_dates %>%
  filter(!is.na(stage_4_date)) %>%
  transmute(patid, date = stage_4_date) %>%
  union_all(
    stage45_dates %>% filter(!is.na(stage_5_date)) %>%
      transmute(patid, date = stage_5_date)
  ) %>%
  union_all(
    ckd5_code_dates %>% transmute(patid, date = ckd5_code_date)
  )

# For each patient, choose the earliest of ALL qualifying routes.
# left_join adds the three supporting dates; dates not present remain NULL/NA.
# The *_at_index flags identify the route(s) that supplied the index date.
# Tied dates can give more than one TRUE flag; no arbitrary hierarchy is imposed.
rk_advanced_ckd_ids <- entry_dates %>%
  group_by(patid) %>%
  summarise(index_date = min(date, na.rm = TRUE), .groups = "drop") %>%
  left_join(stage45_dates, by = "patid") %>%
  left_join(ckd5_code_dates, by = "patid") %>%
  mutate(
    stage_4_at_index = coalesce(stage_4_date == index_date, FALSE),
    stage_5_at_index = coalesce(stage_5_date == index_date, FALSE),
    ckd5_code_at_index = coalesce(ckd5_code_date == index_date, FALSE)
  ) %>%
  anti_join(practice_exclusion_ids, by = "patid") %>%
  anti_join(gender_exclusion_ids, by = "patid") %>%
  output_analysis$cached(
    cache_name("advanced_ckd_ids"),
    unique_indexes = "patid", indexes = "index_date"
  )

# 6. Restrict to the diabetes cohort --------------------------------------------
# Retain all columns in diabetes_cohort and add the kidney-eligibility columns.
# The unique patid index checks that the final table has one row per patient.
rk_diabetes_advanced_ckd_cohort <- diabetes_cohort %>%
  inner_join(rk_advanced_ckd_ids, by = "patid") %>%
  output_analysis$cached(
    cache_name("diabetes_advanced_ckd_cohort"),
    unique_indexes = "patid", indexes = "index_date"
  )

# 7. Small summaries: no need to download the full cohort -----------------------
rk_diabetes_advanced_ckd_cohort %>% count() %>% print()

rk_diabetes_advanced_ckd_cohort %>%
  summarise(
    earliest_index = min(index_date, na.rm = TRUE),
    latest_index = max(index_date, na.rm = TRUE)
  ) %>% print(width = Inf)

# Counts by index-date evidence, including tied evidence routes.
rk_diabetes_advanced_ckd_cohort %>%
  count(stage_4_at_index, stage_5_at_index, ckd5_code_at_index) %>%
  print(n = Inf, width = Inf)

# In a LATER script, reload the final cohort without rebuilding this pipeline:
# analysis <- cprd$analysis("ckd")
# rk_diabetes_advanced_ckd_cohort <- rk_diabetes_advanced_ckd_cohort %>%
#   analysis$cached("rk_elig_v1_diabetes_advanced_ckd_cohort")
# Use the corresponding v2 name if you change cache_version.
#
# Completed persistent caches survive R sessions. Caching does not keep R alive
# or guarantee that an interrupted SQL write completed successfully.
