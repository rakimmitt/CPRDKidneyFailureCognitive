library(aurum)
library(tidyverse)

cprd = CPRDData$new(cprdEnv = "diabetes-2020",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")

codesets = cprd$codesets()
codesets$listCodeSets() %>% print(n=300)

#remove old codelist in order to upload new version
codesets$deleteCodeSet("icd10_ckd5_code") 
codesets$deleteCodeSet("icd10_kf_death") 
codesets$deleteCodeSet("ckd5_code") 

#load medcodelists
icd10_ckd5_code = readr::read_tsv(
  here::here("C:/Users/tj358/OneDrive - University of Exeter/CPRD/Aurum codelists/medcodes/exeter_icd10_ckd5_code.txt"),
  col_types = cols(.default=col_character()))

icd10_ckd5_code %>% codesets$loadICD10CodeSet(name = "icd10_ckd5_code",version="31/10/2021",colname="ICD10")

icd10_kf_death = readr::read_tsv(
  here::here("C:/Users/tj358/OneDrive - University of Exeter/CPRD/Aurum codelists/medcodes/exeter_icd10_kf_death.txt"),
  col_types = cols(.default=col_character()))

icd10_kf_death %>% codesets$loadICD10CodeSet(name = "icd10_kf_death",version="31/10/2021", colname="ICD10")

 ckd5_code = readr::read_tsv(
   here::here("C:/Users/tj358/OneDrive - University of Exeter/CPRD/Aurum codelists/medcodes/exeter_medcodelist_ckd5_code.txt"),
   col_types = cols(.default=col_character()))

 ckd5_code %>% codesets$loadMedCodeSet(name = "ckd5_code",version="31/10/2021",colname="medcodeid")

 codesets$listCodeSets() %>% print(n=300)
 