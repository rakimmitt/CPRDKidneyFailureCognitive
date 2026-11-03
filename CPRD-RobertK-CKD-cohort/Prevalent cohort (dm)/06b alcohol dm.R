# Extracts dates for alcohol consumption code occurrences in GP records

# Merges with index dates (diagnosis dates)

# Defines alcohol status at index date according to our algorithm, described here: https://github.com/Exeter-Diabetes/CPRD-Codelists#alcohol-consumption


############################################################################################

# Setup
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix = "rk_ckd"
draft_suffix <- "pen1_eGFR60+"

analysis = cprd$analysis(analysis_prefix)

############################################################################################

# Pull out all raw code instances and cache with 'all_patid' prefix

analysis = cprd$analysis("all_patid")

raw_alcohol_medcodes <- cprd$tables$observation %>%
  inner_join(codes$alcohol, by="medcodeid") %>%
  analysis$cached("raw_alcohol_medcodes", indexes=c("patid", "obsdate"))


############################################################################################

# Clean: remove if before DOB or after lcd/deregistration, and re-cache

clean_alcohol_medcodes <- raw_alcohol_medcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
  select(patid, date=obsdate, alcohol_cat) %>%
  distinct() %>%
  analysis$cached("clean_alcohol_medcodes", indexes=c("patid", "date", "alcohol_cat"))


############################################################################################

# Find alcohol status according to algorithm at index dates

## Get index dates (diagnosis dates)

analysis = cprd$analysis(analysis_prefix)

matched_cohort <- matched_cohort %>% analysis$cached(paste0("matched_cohort", draft_suffix))
cohort_index_dates <- matched_cohort %>% distinct(patid, index_date)

## Join with alcohol codes on patid and retain codes before index date or up to 7 days after

analysis = cprd$analysis("rk_ckd")

pre_index_date_alcohol_codes <- cohort_index_dates %>%
  inner_join(clean_alcohol_medcodes, by="patid") %>%
  filter(datediff(date, index_date)<=7) %>%
  analysis$cached("alc_pre", indexes = list(c("patid", "index_date")))

## Find if ever previously a 'harmful' drinker (category 3)

harmful_drinker_ever <- pre_index_date_alcohol_codes %>%
  filter(alcohol_cat=="AlcoholConsumptionLevel3") %>%
  distinct(patid, index_date) %>%
  mutate(harmful_drinker_ever=1L)

## Find most recent code
### If different categories on same day, use highest

most_recent_code <- pre_index_date_alcohol_codes %>%
  mutate(alcohol_cat_numeric = case_when(
      alcohol_cat == "AlcoholConsumptionLevel0" ~ 0L,
      alcohol_cat == "AlcoholConsumptionLevel1" ~ 1L,
      alcohol_cat == "AlcoholConsumptionLevel2" ~ 2L,
      alcohol_cat == "AlcoholConsumptionLevel3" ~ 3L,
      TRUE ~ NA_integer_
    )
  ) %>%
  group_by(patid, index_date) %>%
  filter(date==max(date, na.rm=TRUE)) %>%
  summarise(alcohol_cat_numeric = max(alcohol_cat_numeric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  analysis$cached("alc_recent", unique_indexes = list(c("patid", "index_date"))
  )

## Pull together

alcohol_cat <- cohort_index_dates %>%
  left_join(harmful_drinker_ever, c("patid", "index_date")) %>%
  left_join(most_recent_code, c("patid", "index_date")) %>%
  mutate(alcohol_cat_numeric=ifelse(!is.na(harmful_drinker_ever) & harmful_drinker_ever==1, 3L, alcohol_cat_numeric),
         alcohol_cat=case_when(
           alcohol_cat_numeric==0 ~ "None",
           alcohol_cat_numeric==1 ~ "Within limits",
           alcohol_cat_numeric==2 ~ "Excess",
           alcohol_cat_numeric==3 ~ "Harmful",
           TRUE ~ NA_character_
         )) %>%
  select(patid, index_date, alcohol_cat) %>%
  analysis$cached("alcohol", unique_indexes=list(c("patid", "index_date")))
  