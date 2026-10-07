############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix <- "rk_ckd"

############################################################################################

## Cohort and patient characteristics
analysis = cprd$analysis("all")
death_causes <- death_causes %>% analysis$cached("death_causes")

analysis = cprd$analysis("all_patid")
townsend_score <- townsend_score %>% analysis$cached("townsend_score")

analysis = cprd$analysis(analysis_prefix)
ckd_causes <- ckd_causes %>% analysis$cached("ckd_causes")

matched_cohort <- matched_cohort %>% analysis$cached("matched_cohort", unique_indexes="patid")

# create empty dataframe for counts of total population / subset with CKD
counts <- data.frame()

## Biomarkers plus CKD stage
ckd_stages <- ckd_stages %>% analysis$cached("ckd_stages") #this may need to be amended
baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")
  
## Comorbidities
comorbidities <- comorbidities %>% analysis$cached("comorbidities")
  
## Smoking status
smoking <- smoking %>% analysis$cached("smoking")
  
## Medications
medications <- medications %>% analysis$cached("medications")

# EFI

efi <- efi %>% analysis$cached("efi")
  
  
############################################################################################

# Final merge of all datasets to create final cohort for analysis
  
  rk_final_merge <- matched_cohort %>%
    left_join(ckd_stages, by="patid") %>%
    left_join(baseline_biomarkers, by="patid") %>%
    left_join(comorbidities, by="patid") %>%
    left_join(ckd_causes, by="patid") %>%
    left_join(smoking, by="patid") %>%
    left_join(medications, by="patid") %>%
    left_join(townsend_score %>% select(patid, tds_2011), by = "patid") %>% 
    left_join(efi %>% select(patid, efi_n_deficits, pre_index_date_efi_score, pre_index_date_efi_cat), by = "patid") %>%
    left_join(death_causes, by = "patid") %>%
    mutate(index_date_age=datediff(index_date, dob)/365.25,
           index_date_ckd_dur_all=datediff(index_date, first_ckd_date)/365.25,
           index_date = index_date) %>%
    relocate(c(index_date_age, index_date_ckd_dur_all), .before=gender) %>%
    analysis$cached("rk_final_merge", unique_indexes="patid")
  
  ############################################################################################
  
# Export to R data object
# Preserve all integer64 identifiers exactly as character strings.

prev_cohort <- rk_final_merge %>%
  collect() %>%
  mutate(
    across(where(bit64::is.integer64), as.character),
    index_date = as.Date(index_date)
  )

today <- format(Sys.Date(), "%d%m%Y")

save(prev_cohort, file = paste0("C:/Users/rk535/OneDrive - University of Exeter/","CPRD/2024/Raw data/", today, "_matched_ckd_cohort_nondm.Rda"))
  
  rm(medications)
  rm(baseline_biomarkers)
  rm(comorbidities)
  rm(ckd_stages)
  rm(smoking)
  rm(rk_final_merge)
