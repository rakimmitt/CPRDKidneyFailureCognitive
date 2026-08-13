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

#icd10 codes

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

icd10_alzheimers_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_alzheimers_code.txt"),
  col_types = cols(.default=col_character()))

icd10_alzheimers_code %>% codesets$loadICD10CodeSet(name = "icd10_alzheimers_code",version="31/10/2021",colname="ICD10")

icd10_vasculardementia_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_vasculardementia_code.txt"),
  col_types = cols(.default=col_character()))

icd10_vasculardementia_code %>% codesets$loadICD10CodeSet(name = "icd10_vasculardementia_code",version="31/10/2021",colname="ICD10")

icd10_otherdementia_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_otherdementia_code.txt"),
  col_types = cols(.default=col_character()))

icd10_otherdementia_code %>% codesets$loadICD10CodeSet(name = "icd10_otherdementia_code",version="31/10/2021",colname="ICD10")

icd10_alldementia_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_alldementia_code.txt"),
  col_types = cols(.default=col_character()))

icd10_alldementia_code %>% codesets$loadICD10CodeSet(name = "icd10_alldementia_code",version="31/10/2021",colname="ICD10")

icd10_mci_code = readr::read_tsv(
  here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_icd10_mci_code.txt"),
  col_types = cols(.default=col_character()))

icd10_mci_code %>% codesets$loadICD10CodeSet(name = "icd10_mci_code",version="31/10/2021",colname="ICD10")

#medcodes

haemodialysis_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_haemodialysis.txt"),
   col_types = cols(.default=col_character()))

haemodialysis_code %>% codesets$loadMedCodeSet(name = "haemodialysis_code",version="31/10/2021",colname="medcodeid")

peritoneal_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_peritoneal.txt"),
   col_types = cols(.default=col_character()))

peritoneal_code %>% codesets$loadMedCodeSet(name = "peritoneal_code",version="31/10/2021",colname="medcodeid")

transplant_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_transplant.txt"),
   col_types = cols(.default=col_character()))

transplant_code %>% codesets$loadMedCodeSet(name = "transplant_code",version="31/10/2021",colname="medcodeid")

eskdnos_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_eskdnos.txt"),
   col_types = cols(.default=col_character()))

eskdnos_code %>% codesets$loadMedCodeSet(name = "eskdnos_code",version="31/10/2021",colname="medcodeid")

alldementia_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_alldementia.txt"),
   col_types = cols(.default=col_character()))

alldementia_code %>% codesets$loadMedCodeSet(name = "alldementia_code",version="31/10/2021",colname="medcodeid")

vasculardementia_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_vasculardementia.txt"),
   col_types = cols(.default=col_character()))

vasculardementia_code %>% codesets$loadMedCodeSet(name = "vasculardementia_code",version="31/10/2021",colname="medcodeid")

alzheimers_code = readr::read_tsv(
   here::here("https://github.com/rakimmitt/CPRDKidneyFailureCognitive/blob/main/exeter_medcodelist_alzheimers.txt"),
   col_types = cols(.default=col_character()))

alzheimers_code %>% codesets$loadMedCodeSet(name = "alzheimers_code",version="31/10/2021",colname="medcodeid")

 codesets$listCodeSets() %>% print(n=300)
 
