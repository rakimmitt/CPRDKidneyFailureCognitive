#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(dplyr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("rk_ckd")

###############################################################################################

# load ckd stages based on egfr only

analysis = cprd$analysis("all_patid")

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>% analysis$cached("ckd_stages_from_algorithm")

# load in acr data
clean_acr_medcodes <- clean_acr_medcodes %>%
  analysis$cached("clean_acr_medcodes", indexes=c("patid", "date", "testvalue"))

clean_acr_from_separate_medcodes <- clean_acr_from_separate_medcodes %>%
  analysis$cached("clean_acr_from_separate_medcodes", indexes=c("patid", "date", "testvalue"))

all_acr <- clean_acr_medcodes %>%
  select(patid, date, testvalue) %>%
  union_all(clean_acr_from_separate_medcodes %>%
              select(patid, date, testvalue))

# select those with acr >=3 mg/mmol
acr_high <- all_acr %>%
  filter(testvalue >= 3)

acr_span <- acr_high %>%
  group_by(patid) %>%
  summarise(
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    n_tests = n()
  )

# confirm 2 readings 3 months apart or longer
confirmed_acr3 <- acr_span %>%
  filter(n_tests >= 2 & datediff(max_date, min_date) >= 90) %>%
  mutate(confirmed_acr3_date = min_date) %>%
  select(patid, confirmed_acr3_date) %>%
  analysis$cached("confirmed_acr3", indexes = c("patid", "confirmed_acr3_date"))

# join with ckd stage
ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  left_join(confirmed_acr3, by = "patid") %>%
  mutate(
    stage_1 = case_when(
      is.na(stage_1) ~ sql("NULL"),
      is.na(confirmed_acr3_date) ~ sql("NULL"),
      confirmed_acr3_date <= stage_1 ~ stage_1,
      confirmed_acr3_date > stage_1  ~ confirmed_acr3_date
    ),
    stage_2 = case_when(
      is.na(stage_2) ~ sql("NULL"),
      is.na(confirmed_acr3_date) ~ sql("NULL"),
      confirmed_acr3_date <= stage_2 ~ stage_2,
      confirmed_acr3_date > stage_2  ~ confirmed_acr3_date
    )
  )  %>% 
  filter(!(is.na(stage_1) & is.na(stage_2) & is.na(stage_3a) & 
             is.na(stage_3b) & is.na(stage_4) & is.na(stage_5))) %>%
  analysis$cached("ckd_stages_from_algorithm_with_acr",
                  indexes = c("patid"))

# create list of ids of people with ckd and combine with diabetes cohort

# Load exclusions
analysis <- cprd$analysis("diabetes_cohort")

practice_exclusion_ids <- practice_exclusion_ids %>%
  analysis$cached("practice_exclusion_ids")

gender_exclusion_ids <- gender_exclusion_ids %>%
  analysis$cached("gender_exclusion_ids")


# Load diabetes cohort, including diagnosis date
analysis <- cprd$analysis("all")

diabetes_cohort <- diabetes_cohort %>%
  analysis$cached("diabetes_cohort")


# Derive CKD dates without assigning case/control status
analysis <- cprd$analysis("rk_ckd")

ckd_matching_dates <- ckd_stages_from_algorithm %>%
  transmute(
    patid,

    # Earliest recorded stage 3a or 3b date
    ckd_stage_3_start_date = as.Date(
      case_when(
        is.na(stage_3a) ~ stage_3b,
        is.na(stage_3b) ~ stage_3a,
        TRUE ~ pmin(stage_3a, stage_3b)
      )
    ),

    # Earliest stage 4 or 5 date, including diagnostic CKD5
    advanced_ckd_index_date = as.Date(
      case_when(
        is.na(stage_4) ~ stage_5,
        is.na(stage_5) ~ stage_4,
        TRUE ~ pmin(stage_4, stage_5)
      )
    )
  ) %>%
  analysis$cached(
    "ckd_matching_dates",
    unique_indexes = "patid"
  )



##################################################

# Join with relevant tables to get demographics, dementia Dx and other data for ckd cohort, advanced ckd cohort and non-ckd cohort

analysis <- cprd$analysis("all_patid")

# Dementia Dx

raw_alldementia_medcodes <- raw_alldementia_medcodes %>% analysis$cached("raw_alldementia_medcodes")
raw_alldementia_icd10 <- raw_alldementia_icd10 %>% analysis$cached("raw_alldementia_icd10")

# Combine Medcodes and ICD10 records
earliest_all_dementia <- raw_alldementia_medcodes %>%
  select(patid, date = obsdate) %>%
  mutate(source = "gp") %>%
  union_all(raw_alldementia_icd10 %>% select(patid, date = epistart) %>% mutate(source = "hes")) %>%
  inner_join(cprd$tables$validDateLookup, by = "patid") %>%
  filter(date >= min_dob, (source == "gp" & date <= gp_end_date) | (source == "hes" & date <= as.Date("2023-03-31"))) %>%
  group_by(patid) %>%
  summarise(earliest_all_dementia = min(date, na.rm = TRUE),
    .groups = "drop") %>%
  analysis$cached("earliest_all_dementia", unique_indexes = "patid")

## DOB

analysis = cprd$analysis("all")

dob <- cprd$tables$observation %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob) %>%
  group_by(patid) %>%
  summarise(earliest_medcode=min(obsdate, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("earliest_medcode", unique_indexes="patid")

dob %>% count() # should be around 45 million

## No-one has missing dob or earliest_medcode so pmin (runs as 'LEAST' in MySQL) works

dob <- dob %>%
  inner_join(cprd$tables$patient, by="patid") %>%
  mutate(dob=as.Date(ifelse(is.na(mob), paste0(yob,"-06-30"), paste0(yob, "-",mob,"-15")))) %>%
  inner_join(cprd$tables$validDateLookup, by = "patid") %>%
  mutate(dob=pmin(dob, earliest_medcode, na.rm=TRUE)) %>%
  mutate(dob=ifelse(regstartdate>=min_dob & regstartdate<dob, regstartdate, dob)) %>%
  select(patid, dob = dob, mob, yob, regstartdate) %>%
  analysis$cached("dob", unique_indexes="patid")

dob <- dob %>% analysis$cached("dob", unique_indexes="patid")

analysis = cprd$analysis("all_patid")
ethnicity <- ethnicity %>% analysis$cached("ethnicity", unique_indexes="patid")

# Get list of all ids

all_ids <- dob %>%
  anti_join(practice_exclusion_ids, by="patid") %>% 
  anti_join(gender_exclusion_ids, by="patid") %>%
  left_join((cprd$tables$patient %>% select(patid, gender, regenddate, pracid)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd, region)), by="pracid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, reg_date_of_death)), by="patid") %>%
  left_join((cprd$tables$patientImd %>% select(patid, imd_decile)), by="patid") %>%
  left_join((cprd$tables$validDateLookup %>% select(patid, gp_end_date)), by="patid") %>%
  left_join((cprd$tables$patidsWithLinkage %>% mutate(with_hes=1L) %>% select(patid, with_hes, hes_end_date)), by="patid") %>%
  mutate(with_hes=ifelse(is.na(with_hes), 0L, 1L)) %>%
  left_join(ethnicity, by="patid") %>%
  select(patid, gender, dob, pracid, prac_region=region, ethnicity_5cat, ethnicity_16cat, ethnicity_qrisk2, imd_decile, regstartdate, gp_end_date, death_date=reg_date_of_death, with_hes, hes_end_date) %>%
  analysis$cached("all_ids", unique_indexes="patid", indexes=c("gender", "dob"))

all_ids %>% count() #44,363,638

# Join ids with dob and other data for CKD cohort to produce matching_pool

analysis <- cprd$analysis("rk_ckd")

matching_pool <- diabetes_cohort %>%
  select(patid, dm_diag_date_all) %>%

  # Restrict to patients with demographic data and no
  # practice/gender exclusion, as defined in all_ids
  inner_join(all_ids, by = "patid") %>%

  # Preserve people without recorded CKD or dementia
  left_join(ckd_matching_dates, by = "patid") %>%
  left_join(earliest_all_dementia, by = "patid") %>%

  # Requirements that do not depend on a case's index date
  filter(
    with_hes == 1L,
    !is.na(dm_diag_date_all)
  ) %>%

  select(
    patid,
    gender,
    dob,
    pracid,
    prac_region,
    ethnicity_5cat,
    ethnicity_16cat,
    ethnicity_qrisk2,
    imd_decile,
    regstartdate,
    gp_end_date,
    death_date,
    with_hes,
    hes_end_date,
    dm_diag_date_all,
    earliest_all_dementia,
    ckd_stage_3_start_date,
    advanced_ckd_index_date
  ) %>%

  analysis$cached(
    "matching_pool",
    unique_indexes = "patid",
    indexes = c(
      "pracid",
      "advanced_ckd_index_date",
      "ckd_stage_3_start_date"
    )
  )

matching_pool %>% count()
