############################################################################################

#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes_2024 = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix = "rk_ckd"

############################################################################################

## Cohort and patient characteristics
analysis = cprd$analysis("all")
ckd_cohort <- ckd_cohort %>% analysis$cached("diabetes_ckd_cohort")
diabetes_cohort <- diabetes_cohort %>% analysis$cached("diabetes_cohort")
death_causes <- death_causes %>% analysis$cached("death_causes")

townsend_analysis <- cprd$analysis("all_patid")
townsend_score <- townsend_score %>%
  townsend_analysis$cached("townsend_score")

analysis = cprd$analysis(analysis_prefix)
matched_cohort <- matched_cohort %>% analysis$cached("matched_cohort", unique_indexes="patid")

## Get index date

analysis = cprd$analysis("rk_ckd")

#advanced_ckd_ids <- advanced_ckd_ids %>% analysis$cached("advanced_ckd_ids", unique_indexes="patid")
#advanced_ckd_ids <- advanced_ckd_ids %>% select(patid, index_date)

matched_cohort <- matched_cohort %>% analysis$cached("matched_cohort", unique_indexes="patid")
matched_cohort <- matched_cohort %>% select(patid, index_date)

# create empty dataframe for counts of total population / subset with CKD
counts <- data.frame()

## Biomarkers and CKD stages

ckd_stages <- ckd_stages %>% analysis$cached("ckd_stages") #this may need to be amended
baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")

## Comorbidities
comorbidities <- comorbidities %>% analysis$cached("comorbidities")
  
## Smoking status
smoking <- smoking %>% analysis$cached("smoking")
  
## Medications
medications <- medications %>% analysis$cached("medications")
  
############################
  
# Final merge of all datasets to create final cohort for analysis
  
  final_merge <- matched_cohort %>%
    left_join(ckd_stages, by="patid") %>%
    left_join(baseline_biomarkers, by="patid") %>%
    left_join(comorbidities, by="patid") %>%
    left_join(smoking, by="patid") %>%
    left_join(medications, by="patid") %>%
    left_join(townsend_score %>% select(patid, tds_2011), by = "patid") %>%
    left_join(death_causes, by = "patid") %>%
    mutate(index_date_age=datediff(index_date, dob)/365.25,
           index_date_ckd_dur_all=datediff(index_date, first_ckd_date)/365.25,
           dm_dur_all=datediff(index_date, dm_diag_date_all)/365.25,
           index_date = index_date) %>%
    relocate(c(index_date_age, index_date_ckd_dur_all), .before=gender) %>%
    analysis$cached("rk_final_merge", unique_indexes="patid")
  
  ############################################################################################
  
  # Export to R data object
  ## Convert integer64 datatypes to double
  
  prev_cohort <- collect(rk_final_merge %>% mutate(patid=as.character(patid)))
  
  is.integer64 <- function(x){
    class(x)=="integer64"
  }
  
  prev_cohort <- prev_cohort %>%
    mutate_if(is.integer64, as.integer) %>%
    mutate(index_date = as.Date(d))
  
  # Create a valid name (no dashes)
  df_name <- paste0("prev_", gsub("-", "_", d), "_dm")
  
  # Assign name
  assign(df_name, prev_cohort, envir = .GlobalEnv)
  
  today <- format(Sys.Date(), "%Y%m%d")
  
  setwd("C:/Users/rk535/OneDrive - University of Exeter/CPRD/2024/Raw data/")
  save(list = df_name, file=paste0(today, "_prev_ckd_cohort_dm_", d, ".Rda"))
  
  rm(medications)
  rm(baseline_biomarkers)
  rm(comorbidities)
  rm(ckd_stages)
  rm(smoking)
  rm(cohort_ids)
  rm(final_merge)
