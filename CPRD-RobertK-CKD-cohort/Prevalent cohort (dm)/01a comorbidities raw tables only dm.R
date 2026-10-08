#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix = "rk_ckd"

############################################################################################

comorbids <- c("acutepancreatitis",
               "af",
               "angina",
               "anxiety_disorders",
               "asthma",
               "bph",
               "bronchiectasis",
               "chronicpancreatitis",
               "ckd5_code",
               "cld",
               "copd",
               "cysticfibrosis",
               "dementia",
               "diabeticnephropathy",
               "dka",
               "falls",
               "fh_premature_cvd", #family history of premature CVD - for QRISK2
               "frailty_simple",
               "haem_cancer",
               "heartfailure",
               "hosp_cause_majoramputation",
               "hosp_cause_minoramputation",
               "hypertension",
               "ihd", #ischaemic heart disease
               "incident_mi",
               "incident_stroke",
               "lowerlimbfracture",
               "micturition_control",
               "myocardialinfarction",
               #"neuropathy",
               "osteoporosis",
               "otherneuroconditions",
               "pad", #peripheral arterial disease
               "pulmonaryfibrosis",
               "pulmonaryhypertension",
               "severe_retinopathy",
               "non_severe_retinopathy",
               "revasc", #revascularisation procedure
               "rheumatoidarthritis",
               "solid_cancer",
               "solidorgantransplant",
               "stroke",
               "tia",  #transient ischaemic attack
               "ukpds_photocoagulation",
               "unstableangina",
               "urinary_frequency",
               #"vitreoushemorrhage",
               "volume_depletion",
               "genital_infection",
               "genital_infection_nonspec",
               "alldementia",
               "alzheimers",
               "ckd5_nokrt",
               "delirium",
               "haemodialysis",
               "mci",
               "peritoneal_dialysis",
               "renalaccessinfection",
               "transplant",
               "vascular_dementia",
               "uti",
               "skininfection",
               "respiratorytractinfection"
)

############################################################################################

analysis = cprd$analysis("all_patid")

for (i in comorbids) {
  
  if (!i %in% c("alldementia",
               "alzheimers",
               "ckd5_nokrt",
               "delirium",
               "haemodialysis",
               "mci",
               "peritoneal_dialysis",
               "renalaccessinfection",
               "transplant",
               "vascular_dementia",
               "uti",
               "skininfection",
               "respiratorytractinfection")) {
     
    if (length(codes[[i]]) > 0) {
      print(paste("making", i, "medcode table"))
      
      raw_tablename <- paste0("raw_", i, "_medcodes")
      
      data <- cprd$tables$observation %>%
        inner_join(codes[[i]], by="medcodeid") %>%
        analysis$cached(raw_tablename, indexes=c("patid", "obsdate"))
      
      assign(raw_tablename, data)
      
    }
    
    if (i!="hypertension" && length(codes[[paste0("icd10_", i)]]) > 0) {
      print(paste("making", i, "ICD10 code table"))
      
      raw_tablename <- paste0("raw_", i, "_icd10")
      
      data <- cprd$tables$hesDiagnosisEpi %>%
        inner_join(codes[[paste0("icd10_",i)]], sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>%
        analysis$cached(raw_tablename, indexes=c("patid", "epistart"))
      
      assign(raw_tablename, data)
      
    }
    
    if (length(codes[[paste0("opcs4_", i)]]) > 0) {
      print(paste("making", i, "OPCS4 code table"))
      
      raw_tablename <- paste0("raw_", i, "_opcs4")
      
      data <- cprd$tables$hesProceduresEpi %>%
        inner_join(codes[[paste0("opcs4_",i)]], sql_on="LHS.OPCS LIKE CONCAT(opcs4,'%')") %>%
        analysis$cached(raw_tablename, indexes=c("patid", "evdate"))
      
      assign(raw_tablename, data)
      
    }
    
  } else {
    
codelist_root <- paste0(
  "C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/",
  "Github clone/CPRDKidneyFailureCognitive/CPRD-Codelists"
)

          if (file.exists(file.path(codelist_root, "Medcodes", paste0("exeter_medcodelist_", i, ".txt")))) 
            {
      print(paste("making", i, "medcode table"))
      
      raw_tablename <- paste0("raw_", i, "_medcodes")
      
      #placeholder variable that all other codelists have on server
      empty_variable = paste0(i, "_cat")
      
      data <- cprd$tables$observation %>%
        inner_join( readr::read_tsv(
          here::here(paste0("C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Github clone\\CPRDKidneyFailureCognitive\\CPRD-Codelists\\Medcodes\\exeter_medcodelist_", i, ".txt")),
          col_types = cols(.default=col_character())) %>%
            rename(medcodeid=MedCodeId) %>%
            select(medcodeid) %>%
            mutate(!!sym(empty_variable) := NA), 
          by="medcodeid", copy = T) %>%
        analysis$cached(raw_tablename, indexes=c("patid", "obsdate"))
      
      assign(raw_tablename, data)
    }
    
            if (file.exists(file.path(codelist_root, "ICD10", paste0("exeter_icd10_", i, ".txt")))) 
            {
      print(paste("making", i, "ICD10 code table"))
      
      raw_tablename <- paste0("raw_", i, "_icd10")
      empty_variable = paste0("icd10_", i, "_cat")
      
      data <- cprd$tables$hesDiagnosisEpi %>%
        inner_join(
          readr::read_tsv(
            here::here(paste0("C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Github clone\\CPRDKidneyFailureCognitive\\CPRD-Codelists\\ICD10\\exeter_icd10_", i, ".txt")),
            col_types = cols(.default=col_character())) %>% 
            rename(icd10 = ICD10) %>%
            select(icd10) %>%
            mutate(!!sym(empty_variable) := NA),
          , sql_on="LHS.ICD LIKE CONCAT(icd10,'%')", copy = T) %>%
        analysis$cached(raw_tablename, indexes=c("patid", "epistart"))
      
      assign(raw_tablename, data)
            }

            if (file.exists(file.path(codelist_root, "OPCS4", paste0("exeter_opcs4_", i, ".txt")))) 
            {
      print(paste("making", i, "OPCS4 code table"))
      
      raw_tablename <- paste0("raw_", i, "_opcs4")
      empty_variable = paste0("_opcs4", i, "_cat")
      
      data <- cprd$tables$hesProceduresEpi %>%
        inner_join(
          readr::read_tsv(
            here::here(paste0("C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Github clone\\CPRDKidneyFailureCognitive\\CPRD-Codelists\\OPCS4\\exeter_opcs4_", i, ".txt")),
            col_types = cols(.default=col_character())) %>% 
            rename(opcs4 = OPCS4) %>%
            select(opcs4) %>%
            mutate(!!sym(empty_variable) := NA),
          , sql_on="LHS.OPCS LIKE CONCAT(opcs4,'%')", copy = T) %>%
        analysis$cached(raw_tablename, indexes=c("patid", "evdate"))
      
      assign(raw_tablename, data)
      
    }
  }
}

## Add to beginning of list so don't have to remake interim tables when add new comorbidity to end of above list
# Make new primary cause hospitalisation for heart failure, incident MI, and incident stroke comorbidities

raw_primary_hhf_icd10 <- raw_heartfailure_icd10 %>%
  filter(d_order==1) %>%
  analysis$cached("raw_primary_hhf_icd10", indexes=c("patid", "epistart"))

raw_primary_incident_mi_icd10 <- raw_incident_mi_icd10 %>%
  filter(d_order==1) %>%
  analysis$cached("raw_primary_incident_mi_icd10", indexes=c("patid", "epistart"))

raw_primary_incident_stroke_icd10 <- raw_incident_stroke_icd10 %>%
  filter(d_order==1) %>%
  analysis$cached("raw_primary_incident_stroke_icd10", indexes=c("patid", "epistart"))

############################################################################################