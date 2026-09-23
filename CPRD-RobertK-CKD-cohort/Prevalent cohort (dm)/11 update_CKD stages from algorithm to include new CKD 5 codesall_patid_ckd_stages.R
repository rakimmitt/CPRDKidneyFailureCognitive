# Calculate longitudinal CKD stages using our algorithm (https://github.com/Exeter-Diabetes/CPRD-Codelists#ckd-chronic-kidney-disease-stage)

## Combine with CKD5 medcodes
## Find start date for each CKD stage
## Reshape wide to give 1 row per patid with start dates of each CKD stage


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("all_patid")

################################################

# Combine with CKD5 medcodes/ICD10/OPCS4 codes

## Get raw CKD5 codes and clean
### All are already in all_patid tables on MySQL from other script

### Medcodes
raw_ckd5_code_medcodes <- raw_ckd5_medcodes %>% analysis$cached("raw_ckd5_code_medcodes")

### ICD10 codes
raw_ckd5_code_icd10 <- raw_ckd5_icd10 %>% analysis$cached("raw_ckd5_code_icd10")

### OPCS4 codes
raw_ckd5_code_opcs4 <- raw_ckd5_opcs4 %>% analysis$cached("raw_ckd5_code_opcs4")


## Clean, find earliest date per person, and re-cache

earliest_clean_ckd5 <- raw_ckd5_code_medcodes %>%
  select(patid, date=obsdate) %>%
  mutate(source="gp") %>%
  union_all((raw_ckd5_code_icd10 %>% select(patid, date=epistart) %>% mutate(source="hes"))) %>%
  union_all((raw_ckd5_code_opcs4 %>% select(patid, date=evdate) %>% mutate(source="hes"))) %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(date>=min_dob & ((source=="gp" & date<=gp_end_date) | (source=="hes" & date<=as.Date("2023-03-31")))) %>%
  group_by(patid) %>%
  summarise(first_test_date=min(date, na.rm=TRUE))%>%
  ungroup() %>%
  analysis$cached("earliest_clean_ckd5",indexes=c("patid", "first_test_date"))


## Combine CKD5 and other codes

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  select(patid, ckd_stage, first_test_date) %>%
  union_all(earliest_clean_ckd5 %>% mutate(ckd_stage="stage_5")) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_5",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()        
#6,123,490


################################################################################################################################

# Define date of onset for each stage

## For each person, define date of onset of each stage (earliest incident) - assume no returning to less severe stages

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  group_by(patid, ckd_stage) %>%
  summarise(ckd_stage_start=min(first_test_date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("ckd_stages_from_algorithm_interim_6",indexes=c("patid","ckd_stage","ckd_stage_start"))
  

## Remove where start date of less severe stage is later than start date of more severe stage
### Reshape wide first

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  pivot_wider(id_cols=patid,
              names_from=ckd_stage,
              values_from=ckd_stage_start) %>%
  mutate(stage_1=ifelse(!is.na(stage_1) & !is.na(stage_2) & stage_1>stage_2, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_3a) & stage_1>stage_3a, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_3b) & stage_1>stage_3b, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_4) & stage_1>stage_4, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_5) & stage_1>stage_5, NA, stage_1),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_3a) & stage_2>stage_3a, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_3b) & stage_2>stage_3b, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_4) & stage_2>stage_4, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_5) & stage_2>stage_5, NA, stage_2),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_3b) & stage_3a>stage_3b, NA, stage_3a),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_4) & stage_3a>stage_4, NA, stage_3a),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_5) & stage_3a>stage_5, NA, stage_3a),
         stage_3b=ifelse(!is.na(stage_3b) & !is.na(stage_4) & stage_3b>stage_4, NA, stage_3b),
         stage_3b=ifelse(!is.na(stage_3b) & !is.na(stage_5) & stage_3b>stage_5, NA, stage_3b),
         stage_4=ifelse(!is.na(stage_4) & !is.na(stage_5) & stage_4>stage_5, NA, stage_4)) %>%
  analysis$cached("ckd_stages_from_algorithm", unique_indexes="patid")
                  
