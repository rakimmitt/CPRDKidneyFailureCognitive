library(aurum)
library(tidyverse)

cprd = CPRDData$new(cprdEnv = "diabetes-2020",cprdConf = "C:/Users/rakim/Documents/.aurum.yaml")

codesets = cprd$codesets()
codesets$listCodeSets() %>% print(n=300)

#remove old codelist in order to upload new version
codesets$deleteCodeSet("icd10_ckd5_code") 
codesets$deleteCodeSet("icd10_kf_death") 
codesets$deleteCodeSet("ckd5_code") 

#load medcodelists

icd10_haemodialysis_sensitive_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_haemodialysis_sensitive_code.txt"),
  col_types = cols(.default=col_character()))

icd10_haemodialysis_sensitive_code %>% codesets$loadICD10CodeSet(name = "icd10_haemodialysis_sensitive_code",version="31/10/2021",colname="ICD10")

icd10_peritoneal_sensitive_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_peritoneal_sensitive_code.txt"),
  col_types = cols(.default=col_character()))

icd10_peritoneal_sensitive_code %>% codesets$loadICD10CodeSet(name = "icd10_peritoneal_sensitive_code",version="31/10/2021",colname="ICD10")

icd10_transplant_sensitive_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_transplant_sensitive_code.txt"),
  col_types = cols(.default=col_character()))

icd10_transplant_sensitive_code %>% codesets$loadICD10CodeSet(name = "icd10_transplant_sensitive_code",version="31/10/2021",colname="ICD10")

icd10_eskdnos_sensitive_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_eskdnos_sensitive_code.txt"),
  col_types = cols(.default=col_character()))

icd10_eskdnos_sensitive_code %>% codesets$loadICD10CodeSet(name = "icd10_eskdnos_sensitive_code",version="31/10/2021",colname="ICD10")

haemodialysis_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_haemodialysis.txt"),
   col_types = cols(.default=col_character()))

haemodialysis_code %>% codesets$loadMedCodeSet(name = "haemodialysis_code",version="31/10/2021",colname="medcodeid")

 codesets$listCodeSets() %>% print(n=300)
 
